# Phys 39 Module 6, Part 7 - P and PI temperature control GUI
# (copied from the Module 5 P-control GUI, module5_p_control_gui.py)
#
# Changes from Module 5:
#   - PI control law with an integral state q and conditional-integration
#     anti-windup (pi_control() below); Ki = 0 gives the Module 5 P controller.
#   - dt is the actual elapsed time between accepted temperature readings
#     (Arduino time stamps), not a nominal 1 s.
#   - Zero integral button (q = 0 before a new comparison; only while control
#     is off, i.e. with PWM 0). A run cannot start until q = 0.
#   - Two fixed setpoints (low and high) for the Goldilocks runs; the active
#     setpoint can be switched during a run, so one run file holds the upward
#     and the downward step.
#   - Adjustable PWM limit u_max (at most 255) and anti-windup switch, both
#     fixed for a run.
#   - e, u_P, u_I, the requested and the applied command, q, and dt are shown
#     and saved.
#   - Plots keep the whole session (no rolling window). Mouse wheel zooms,
#     drag pans, and "Show all" returns to the full auto-scaled view.
#   - The serial link, reader thread, and line parser are copied in from the
#     Module 4 GUI (with the Module 5 SerialLink fix), so this file no longer
#     imports module4_control_gui.py (now in archive/).
#
# Use with arduino/module6_pi_control/module6_pi_control.ino (9600 baud). Its
# temperature shutdown, 0-255 clamp, and 3 s command timeout are unchanged.
#
# Arduino emits lines like (one per 1000-reading average, about 1 s):
#   Temperature (C): 27.73, Time (s): 645.06, PWM: 120, Direction input: 1, Active PWM pin: 9, Heat/Cool: 1, Safety shutdown: 0
#
# Two modes:
#   MANUAL    - same as Module 4/5: direction drop-down + PWM slider.
#   CONTROL   - on every new averaged temperature T_n, Python calculates
#                   e_n     = Tset - T_n
#                   q_{n+1} = q_n + e_n dt         (held if anti-windup applies)
#                   u_P     = Kp e_n,   u_I = Ki q_{n+1}
#                   u_req   = u_P + u_I            (requested command)
#                   u       = clamp(u_req, -u_max, +u_max)   (applied command)
#               and sends  SET PWM <P> DIR HEAT  if u >= 0
#                          SET PWM <P> DIR COOL  if u <  0
#               with P = |u| rounded to an integer.
#
# The ONLY control calculation is in pi_control() below. Kp, Ki, u_max, and
# the anti-windup setting are fixed for a whole run. No feed-forward, gain
# scheduling, or derivative term. The clamp and the stop-on-Arduino-safety-
# shutdown rule are safety constraints, not part of the control law.
#
# Each run is saved to its own CSV in data/module_06/:
#   <p|pi>_<YYYYmmdd-HHMMSS>_Kp<gain>_Ki<gain>_T<low>-<high>.csv
# (with _umax<limit> if u_max < 255 and _noAW if anti-windup is off)
# columns:
#   time_s, dt_s, temperature_C, setpoint_C, Kp, Ki, u_max, anti_windup,
#   error_C, q_Cs, u_P, u_I, u_req, u_applied, integral_held, saturated,
#   pwm_cmd, dir_cmd, pwm_arduino, heat_cool, safety_shutdown, control_on
# u_applied/pwm_cmd/dir_cmd are what Python sent after this temperature;
# pwm_arduino/heat_cool are what the Arduino reported it was applying when it
# printed this line (the previous command), so they lag by one line.
#
# As in Module 5, each run file starts with the last BASELINE_SECONDS of
# readings from BEFORE control was switched on (control_on = 0, nothing
# sent by the controller), then every reading while control is on.

import re
import sys
import threading
from collections import deque
from datetime import datetime
from typing import NamedTuple

import serial
from PySide6 import QtCore, QtWidgets
import pyqtgraph as pg

from module6_common import MODULE6_DATA_DIR, TAU_S, module4_susceptibility, one_lump_pi_zeta


# ------------------------------------------------------------------
# USER CONFIGURATION  (edit these near the top)
# ------------------------------------------------------------------
SERIAL_PORT = "COM5"          # "/dev/cu.usbmodem1101" Mac. Windows example: "COM5"
BAUD_RATE = 9600

PLOT_UPDATE_MS = 100          # how often the plots redraw (ms)
SEND_DEBOUNCE_MS = 80         # manual mode: wait after the last change before sending
TEMP_MIN_Y_SPAN = 0.2         # smallest temperature-axis span (C): limits how far the
                              # noise is magnified, still lets you zoom to the steady error

PWM_MAX = 255                 # 8-bit PWM magnitude limit (also clamped by the Arduino)

# Allowed setpoint range in the GUI and the two Goldilocks setpoints.
# Use the same two setpoints for all six runs.
SETPOINT_MIN_C = 10.0
SETPOINT_MAX_C = 50.0
SETPOINT_LOW_DEFAULT_C = 25.0
SETPOINT_HIGH_DEFAULT_C = 30.0

KP_MAX = 100.0                # largest Kp the GUI accepts (PWM per C)
KP_DEFAULT = 0.5
KI_MAX = 5.0                  # largest Ki the GUI accepts (PWM per C s)
KI_DEFAULT = 0.0              # P only until Ki is set deliberately

DATA_DIR = MODULE6_DATA_DIR
BASELINE_SECONDS = 30         # readings before Start written to each run file
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# PI CONTROL LAW  (the lines to show the instructor before powering)
# ------------------------------------------------------------------
class PIStep(NamedTuple):
    e: float            # error e_n = Tset - T_n                      [C]
    q: float            # integral state q_{n+1} after this step      [C s]
    u_p: float          # u_P = Kp e_n                                [PWM]
    u_i: float          # u_I = Ki q_{n+1}                            [PWM]
    u_req: float        # requested command u_P + u_I                 [PWM]
    u_applied: float    # clamped command, -u_max ... +u_max          [PWM]
    held: bool          # True if conditional integration kept q_n
    saturated: bool     # True if |u_req| > u_max
    direction: str      # "HEAT" if u_applied >= 0, else "COOL"
    pwm: int            # |u_applied| rounded, sent to the Arduino     [0-255]


def pi_control(t_set: float, temperature: float, q: float, dt: float,
               kp: float, ki: float, u_max: float, anti_windup: bool) -> PIStep:
    """One PI update for one accepted temperature (Module 6 handout, Part 7):
        q_{n+1} = q_n + e_n dt,    u_n = Kp e_n + Ki q_{n+1}.
    Ki = 0 gives P-only control (u_I = 0)."""
    e = t_set - temperature                      # error e_n = Tset - T_n           [C]

    # Integral state. Conditional integration (anti-windup): if the request
    # with the updated q is beyond the PWM limit and this error pushes it
    # farther, keep q_n on this step. Otherwise update q, including when the
    # error helps bring the command back out of saturation.
    q_next = q + e * dt                          # q_{n+1} = q_n + e_n dt           [C s]
    u_next = kp * e + ki * q_next
    held = anti_windup and ((u_next > u_max and e > 0) or (u_next < -u_max and e < 0))
    if held:
        q_next = q

    u_p = kp * e                                 # u_P = Kp e_n                     [PWM]
    u_i = ki * q_next                            # u_I = Ki q_{n+1}                 [PWM]
    u_req = u_p + u_i                            # requested command                [PWM]
    u_applied = max(-u_max, min(u_max, u_req))   # clamped (applied) command        [PWM]

    # Convert u to the Arduino command: the sign selects the direction,
    # the magnitude |u| becomes the integer PWM value.
    direction = "HEAT" if u_applied >= 0 else "COOL"
    pwm = min(int(abs(u_applied) + 0.5), PWM_MAX)
    saturated = abs(u_req) > u_max
    return PIStep(e, q_next, u_p, u_i, u_req, u_applied, held, saturated, direction, pwm)


# ------------------------------------------------------------------
# SERIAL LINK  (from the Module 4 GUI, with the Module 5 fix)
# ------------------------------------------------------------------
# The reader thread reads lines and the GUI thread writes commands.
# pyserial allows one thread to read while another writes, so only the
# writes are serialised with a lock. (Module 5 fix: the first version also
# held the lock while the reader waited for a line, which delayed commands.)
class SerialLink:
    def __init__(self, port: str, baud: int):
        self.ser = serial.Serial(port, baud, timeout=0.5)
        self.write_lock = threading.Lock()

    def readline(self) -> bytes:
        # Only the reader thread calls this, so no lock.
        return self.ser.readline()

    def write_line(self, text: str):
        with self.write_lock:
            self.ser.write((text + "\n").encode("ascii"))
            self.ser.flush()

    def close(self):
        # Callers stop the reader thread first (SerialReader.stop()).
        with self.write_lock:
            try:
                self.ser.close()
            except Exception:
                pass


# Parser for one Arduino measurement line. Direction input and Active PWM pin
# are matched but not used; Safety shutdown is optional.
MEAS_RE = re.compile(
    r"Temperature \(C\):\s*(?P<temp>-?\d+(?:\.\d+)?|nan)\s*,\s*"
    r"Time \(s\):\s*(?P<time>-?\d+(?:\.\d+)?)\s*,\s*"
    r"PWM:\s*(?P<pwm>\d+)\s*,\s*"
    r"Direction input:\s*(?P<dir>[01])\s*,\s*"
    r"Active PWM pin:\s*(?P<pin>\d+)\s*,\s*"
    r"Heat/Cool:\s*(?P<hc>[01])"
    r"(?:\s*,\s*Safety shutdown:\s*(?P<safety>[01]))?"
)


def parse_measurement(line: str):
    """Return a dict of the fields, or None if the line is not an accepted reading."""
    m = MEAS_RE.search(line)
    if not m:
        return None
    if m.group("temp").lower() == "nan":
        return None
    return {
        "time_s": float(m.group("time")),
        "temperature_C": float(m.group("temp")),
        "pwm": int(m.group("pwm")),
        "heat_cool": int(m.group("hc")),
        "safety_shutdown": int(m.group("safety") or 0),
    }


class SerialReader(QtCore.QThread):
    line_received = QtCore.Signal(str)
    error = QtCore.Signal(str)

    def __init__(self, link):
        super().__init__()
        self.link = link
        self._stop = False

    def run(self):
        try:
            while not self._stop:
                raw = self.link.readline()
                if not raw:
                    continue
                line = raw.decode("utf-8", errors="replace").strip()
                if line:
                    self.line_received.emit(line)
        except Exception as e:
            self.error.emit(str(e))

    def stop(self):
        self._stop = True
        self.wait(1000)


# ------------------------------------------------------------------
# GUI + PLOTS
# ------------------------------------------------------------------
class ControlGUI(QtWidgets.QMainWindow):
    def __init__(self, link=None):
        super().__init__()
        self.setWindowTitle("Module 6 Part 7 - P/PI TEC Temperature Control")

        # ------- Module 4 susceptibilities (for the gain display) -------
        try:
            self.chi_h, self.chi_c = module4_susceptibility()
        except Exception as e:
            print(f"[warning] could not load Module 4 susceptibility: {e}")
            self.chi_h = self.chi_c = None

        # ------- controller state -------
        self.control_on = False
        self.q = 0.0                 # integral state q (C s); zeroed only by the button
        self.t_prev = None           # Arduino time of the previous accepted reading (s)
        self.t_offset = 0.0          # plot time = Arduino time + offset (grows if the Arduino resets)
        self.run = None              # gains and limits fixed for the current run
        self.csv_file = None
        self.csv_path = None
        # Recent readings while control is off (about one per second), so
        # each run file starts with a baseline before the step.
        self.baseline = deque(maxlen=BASELINE_SECONDS)

        # ------- in-memory data for the plots (whole session, never trimmed) -------
        self.ts = []                 # plot time (s): Arduino time, continued across resets
        self.temps = []              # temperature (C)
        self.setpoints = []          # setpoint (C), NaN when control is off
        self.errors = []             # error (C), NaN when control is off
        self.u_ps = []               # u_P (PWM), NaN when control is off
        self.u_is = []               # u_I (PWM), NaN when control is off
        self.u_reqs = []             # requested command (PWM), NaN when control is off
        self.pwm_signed = []         # Arduino-reported signed PWM (cooling negative)

        # ------- open serial link -------
        self.link = link if link is not None else SerialLink(SERIAL_PORT, BAUD_RATE)

        # ------- build the UI -------
        self._build_ui()

        # ------- start serial reader thread -------
        self.reader = SerialReader(self.link)
        self.reader.line_received.connect(self.on_line)
        self.reader.error.connect(self.on_serial_error)
        self.reader.start()

        # ------- plot update timer -------
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.update_plots)
        self.timer.start(PLOT_UPDATE_MS)

        # ------- manual command send debounce (as in Module 4) -------
        self.send_timer = QtCore.QTimer(self)
        self.send_timer.setSingleShot(True)
        self.send_timer.setInterval(SEND_DEBOUNCE_MS)
        self.send_timer.timeout.connect(self.send_manual_command)

        # Safety: PWM starts at zero (SET PWM 0 DIR HEAT at startup).
        self.schedule_send()

    # ---------------- build the GUI ----------------
    def _build_ui(self):
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QVBoxLayout(central)

        # --- manual controls (Module 4) ---
        manual_box = QtWidgets.QGroupBox("Manual (open loop)")
        row1 = QtWidgets.QHBoxLayout(manual_box)

        row1.addWidget(QtWidgets.QLabel("Direction:"))
        self.dir_combo = QtWidgets.QComboBox()
        self.dir_combo.addItems(["HEAT", "COOL"])
        self.dir_combo.currentIndexChanged.connect(self.schedule_send)
        row1.addWidget(self.dir_combo)

        row1.addSpacing(20)
        row1.addWidget(QtWidgets.QLabel("PWM:"))
        self.pwm_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.pwm_slider.setMinimum(0)
        self.pwm_slider.setMaximum(PWM_MAX)
        self.pwm_slider.setValue(0)
        self.pwm_slider.valueChanged.connect(self.on_slider_changed)
        row1.addWidget(self.pwm_slider, stretch=1)

        self.pwm_edit = QtWidgets.QLineEdit("0")
        self.pwm_edit.setFixedWidth(60)
        self.pwm_edit.editingFinished.connect(self.on_text_edited)
        row1.addWidget(self.pwm_edit)

        layout.addWidget(manual_box)

        # --- P/PI control ---
        pi_box = QtWidgets.QGroupBox("P/PI control:  u = Kp e + Ki q,   e = Tset - T,   q = sum of e dt")
        grid = QtWidgets.QGridLayout(pi_box)

        # setpoints: two fixed values, the active one can change during a run
        grid.addWidget(QtWidgets.QLabel("T low (C):"), 0, 0)
        self.low_spin = self._spin(SETPOINT_MIN_C, SETPOINT_MAX_C, 1, 0.5, SETPOINT_LOW_DEFAULT_C)
        grid.addWidget(self.low_spin, 0, 1)
        grid.addWidget(QtWidgets.QLabel("T high (C):"), 0, 2)
        self.high_spin = self._spin(SETPOINT_MIN_C, SETPOINT_MAX_C, 1, 0.5, SETPOINT_HIGH_DEFAULT_C)
        grid.addWidget(self.high_spin, 0, 3)
        grid.addWidget(QtWidgets.QLabel("Active Tset:"), 0, 4)
        self.low_button = QtWidgets.QPushButton("Tset = low")
        self.high_button = QtWidgets.QPushButton("Tset = high")
        self.setpoint_group = QtWidgets.QButtonGroup(self)
        for b in (self.low_button, self.high_button):
            b.setCheckable(True)
            self.setpoint_group.addButton(b)
        self.high_button.setChecked(True)
        self.setpoint_group.buttonClicked.connect(self.on_setpoint_selected)
        grid.addWidget(self.low_button, 0, 5)
        grid.addWidget(self.high_button, 0, 6)

        # gains, PWM limit, anti-windup
        grid.addWidget(QtWidgets.QLabel("Kp (PWM/C):"), 1, 0)
        self.kp_spin = self._spin(0.0, KP_MAX, 3, 0.1, KP_DEFAULT)
        grid.addWidget(self.kp_spin, 1, 1)
        grid.addWidget(QtWidgets.QLabel("Ki (PWM/(C s)):"), 1, 2)
        self.ki_spin = self._spin(0.0, KI_MAX, 3, 0.01, KI_DEFAULT)
        grid.addWidget(self.ki_spin, 1, 3)
        grid.addWidget(QtWidgets.QLabel("u_max (PWM):"), 1, 4)
        self.umax_spin = QtWidgets.QSpinBox()
        self.umax_spin.setRange(1, PWM_MAX)
        self.umax_spin.setValue(PWM_MAX)
        grid.addWidget(self.umax_spin, 1, 5)
        self.antiwindup_check = QtWidgets.QCheckBox("Anti-windup (conditional integration)")
        self.antiwindup_check.setChecked(True)
        grid.addWidget(self.antiwindup_check, 1, 6)
        for w in (self.kp_spin, self.ki_spin):
            w.valueChanged.connect(self.update_gain_label)

        # one-lump predictions, integral state, buttons
        self.lbl_gain = QtWidgets.QLabel()
        grid.addWidget(self.lbl_gain, 2, 0, 1, 5)
        self.zero_button = QtWidgets.QPushButton("Zero integral (q = 0)")
        self.zero_button.clicked.connect(self.zero_integral)
        grid.addWidget(self.zero_button, 2, 5)
        self.control_button = QtWidgets.QPushButton("Start control")
        self.control_button.clicked.connect(self.on_control_button)
        grid.addWidget(self.control_button, 2, 6)
        self.update_gain_label()

        layout.addWidget(pi_box)

        # --- live value labels ---
        row3 = QtWidgets.QHBoxLayout()
        self.lbl_temp = QtWidgets.QLabel("T: --- C")
        self.lbl_tset = QtWidgets.QLabel("Tset: ---")
        self.lbl_error = QtWidgets.QLabel("e: ---")
        self.lbl_up = QtWidgets.QLabel("u_P: ---")
        self.lbl_ui = QtWidgets.QLabel("u_I: ---")
        self.lbl_ureq = QtWidgets.QLabel("u_req: ---")
        self.lbl_u = QtWidgets.QLabel("applied: ---")
        self.lbl_q = QtWidgets.QLabel("q: 0.0 C s")
        self.lbl_dt = QtWidgets.QLabel("dt: ---")
        self.lbl_state = QtWidgets.QLabel("")
        for lbl in (self.lbl_temp, self.lbl_tset, self.lbl_error, self.lbl_up, self.lbl_ui,
                    self.lbl_ureq, self.lbl_u, self.lbl_q, self.lbl_dt, self.lbl_state):
            row3.addWidget(lbl)
        row3.addStretch(1)
        layout.addLayout(row3)

        row4 = QtWidgets.QHBoxLayout()
        self.lbl_pwm = QtWidgets.QLabel("Arduino PWM: ---")
        self.lbl_dir = QtWidgets.QLabel("Direction: ---")
        self.lbl_time = QtWidgets.QLabel("Time: --- s")
        self.lbl_safety = QtWidgets.QLabel("Safety: ---")
        for lbl in (self.lbl_pwm, self.lbl_dir, self.lbl_time, self.lbl_safety):
            row4.addWidget(lbl)
        row4.addStretch(1)
        row4.addWidget(QtWidgets.QLabel(
            "Plots: wheel = zoom (over an axis: that axis only), drag = pan, right-click = menu"))
        self.show_all_button = QtWidgets.QPushButton("Show all")
        self.show_all_button.clicked.connect(self.show_all)
        row4.addWidget(self.show_all_button)
        layout.addLayout(row4)

        # --- temperature plot with setpoint ---
        self.temp_plot = pg.PlotWidget()
        self.temp_plot.setLabel("left", "Temperature (C)")
        self.temp_plot.addLegend(offset=(10, 10))
        self.temp_curve = self.temp_plot.plot(pen=pg.mkPen("k", width=2), name="T")
        self.setpoint_curve = self.temp_plot.plot(
            pen=pg.mkPen("g", width=2, style=QtCore.Qt.DashLine),
            connect="finite", name="Tset",
        )
        self.temp_plot.getViewBox().setLimits(minYRange=TEMP_MIN_Y_SPAN)
        layout.addWidget(self.temp_plot, stretch=2)

        # --- command plot: u_P, u_I, requested, and the PWM the Arduino applied ---
        self.cmd_plot = pg.PlotWidget()
        self.cmd_plot.setLabel("left", "Signed PWM")
        self.cmd_plot.addLegend(offset=(10, 10))
        # Drawn first, so u_I stays visible where it equals the applied command.
        self.pwm_curve = self.cmd_plot.plot(pen=pg.mkPen("k", width=2),
                                            name="applied (Arduino; cooling negative)")
        self.up_curve = self.cmd_plot.plot(pen=pg.mkPen("b", width=1.5),
                                           connect="finite", name="u_P")
        self.ui_curve = self.cmd_plot.plot(pen=pg.mkPen((255, 140, 0), width=2),
                                           connect="finite", name="u_I")
        self.ureq_curve = self.cmd_plot.plot(
            pen=pg.mkPen((128, 128, 128), width=1.5, style=QtCore.Qt.DashLine),
            connect="finite", name="requested u_P + u_I")
        self.cmd_plot.addLine(y=0, pen=pg.mkPen((160, 160, 160), width=1)).setZValue(-10)
        self.umax_lines = []
        for sign in (1, -1):
            line = pg.InfiniteLine(pos=sign * PWM_MAX, angle=0,
                                   pen=pg.mkPen("r", width=1, style=QtCore.Qt.DotLine))
            self.cmd_plot.addItem(line, ignoreBounds=True)   # keep it out of auto-range
            self.umax_lines.append(line)
        layout.addWidget(self.cmd_plot, stretch=2)

        # --- error plot ---
        self.error_plot = pg.PlotWidget()
        self.error_plot.setLabel("bottom", "Time (s)")
        self.error_plot.setLabel("left", "e = Tset - T (C)")
        self.error_plot.addLine(y=0, pen=pg.mkPen("k", width=1)).setZValue(-10)   # behind e
        self.error_curve = self.error_plot.plot(
            pen=pg.mkPen("m", width=2), connect="finite"
        )
        layout.addWidget(self.error_plot, stretch=1)

        # All three plots share the time axis. Auto-range shows the whole session;
        # zooming or panning one plot turns auto-range off until "Show all".
        self.plots = (self.temp_plot, self.cmd_plot, self.error_plot)
        for p in self.plots:
            p.showGrid(x=True, y=True, alpha=0.3)
            if p is not self.temp_plot:
                p.setXLink(self.temp_plot)

        self.resize(1250, 950)

    def _spin(self, low, high, decimals, step, value):
        spin = QtWidgets.QDoubleSpinBox()
        spin.setRange(low, high)
        spin.setDecimals(decimals)
        spin.setSingleStep(step)
        spin.setValue(value)
        return spin

    def update_gain_label(self):
        kp, ki = self.kp_spin.value(), self.ki_spin.value()
        if self.chi_h is None:
            self.lbl_gain.setText("Module 4 data not found")
            return
        text = (f"L_heat = Kp chi_h = {kp * self.chi_h:.2f},  "
                f"L_cool = Kp chi_c = {kp * self.chi_c:.2f}.   ")
        if ki > 0:
            text += (f"One-lump zeta: heat {one_lump_pi_zeta(self.chi_h, kp, ki):.2f}, "
                     f"cool {one_lump_pi_zeta(self.chi_c, kp, ki):.2f}.   ")
        else:
            text += "Ki = 0: P only.   "
        text += f"No one-lump overshoot if Ki <= Kp/tau = {kp / TAU_S:.3f}"
        self.lbl_gain.setText(text)

    # ---------------- manual control callbacks ----------------
    def on_slider_changed(self, value: int):
        if self.pwm_edit.text() != str(value):
            self.pwm_edit.setText(str(value))
        self.schedule_send()

    def on_text_edited(self):
        text = self.pwm_edit.text().strip()
        try:
            value = int(text)
        except ValueError:
            value = self.pwm_slider.value()
        value = max(0, min(PWM_MAX, value))
        self.pwm_edit.setText(str(value))
        if self.pwm_slider.value() != value:
            self.pwm_slider.setValue(value)   # triggers on_slider_changed
        else:
            self.schedule_send()

    def schedule_send(self):
        self.send_timer.start()

    def send_manual_command(self, echo: bool = True):
        # Manual commands are ignored while the controller owns the output.
        if self.control_on:
            return
        pwm = self.pwm_slider.value()
        direction = self.dir_combo.currentText()
        self.send_line(f"SET PWM {pwm} DIR {direction}", echo)

    def send_line(self, cmd: str, echo: bool = True):
        try:
            self.link.write_line(cmd)
            if echo:
                print(f"[sent] {cmd}")
        except Exception as e:
            print(f"[send error] {e}")

    # ---------------- setpoint and integral ----------------
    def active_setpoint(self) -> float:
        if self.control_on:
            return self.run["low"] if self.low_button.isChecked() else self.run["high"]
        return self.low_spin.value() if self.low_button.isChecked() else self.high_spin.value()

    def on_setpoint_selected(self, _button):
        if self.control_on:
            t = self.t_prev if self.t_prev is not None else float("nan")
            print(f"[setpoint] t = {t:.2f} s: Tset -> {self.active_setpoint():.2f} C")

    def zero_integral(self):
        # Only reachable while control is off (the button is disabled during a
        # run), i.e. with the controller's PWM already at 0.
        self.q = 0.0
        self.lbl_q.setText("q: 0.0 C s")
        print("[integral] q = 0")

    # ---------------- control start / stop ----------------
    def on_control_button(self):
        if self.control_on:
            self.stop_control("stopped by user")
        else:
            self.start_control()

    def start_control(self):
        if self.q != 0.0:
            QtWidgets.QMessageBox.warning(
                self, "Integral not zero",
                f"q = {self.q:.1f} C s from the previous run.\n"
                "Press 'Zero integral' before starting a new run.")
            return
        low, high = self.low_spin.value(), self.high_spin.value()
        if not low < high:
            QtWidgets.QMessageBox.warning(self, "Setpoints",
                                          "T low must be below T high.")
            return
        self.run = {
            "low": low, "high": high,
            "kp": self.kp_spin.value(), "ki": self.ki_spin.value(),
            "u_max": float(self.umax_spin.value()),
            "anti_windup": self.antiwindup_check.isChecked(),
        }
        kp, ki = self.run["kp"], self.run["ki"]

        DATA_DIR.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        kind = "pi" if ki > 0 else "p"
        extra = ("" if self.run["u_max"] == PWM_MAX else f"_umax{self.run['u_max']:g}") + (
            "" if self.run["anti_windup"] else "_noAW")
        self.csv_path = DATA_DIR / f"{kind}_{stamp}_Kp{kp:g}_Ki{ki:g}_T{low:g}-{high:g}{extra}.csv"
        self.csv_file = open(self.csv_path, "w", newline="")
        self.csv_file.write(
            "time_s,dt_s,temperature_C,setpoint_C,Kp,Ki,u_max,anti_windup,"
            "error_C,q_Cs,u_P,u_I,u_req,u_applied,integral_held,saturated,"
            "pwm_cmd,dir_cmd,pwm_arduino,heat_cool,safety_shutdown,control_on\n"
        )
        t_set = low if self.low_button.isChecked() else high
        for t, dt, T, pwm_ard, hc, safety in self.baseline:
            self.write_row(t, dt, T, t_set, PIStep(t_set - T, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                   False, False, "HEAT", 0),
                           pwm_ard, hc, safety, control_on=0)
        self.baseline.clear()
        self.csv_file.flush()

        self.control_on = True
        self.send_timer.stop()
        # Gains, limit, anti-windup, and the two setpoints stay fixed for the run;
        # only the active setpoint (low/high) can be switched.
        for w in (self.low_spin, self.high_spin, self.kp_spin, self.ki_spin,
                  self.umax_spin, self.antiwindup_check, self.zero_button,
                  self.dir_combo, self.pwm_slider, self.pwm_edit):
            w.setEnabled(False)
        for line, sign in zip(self.umax_lines, (1, -1)):
            line.setValue(sign * self.run["u_max"])
        self.control_button.setText("Stop control (PWM 0)")
        print(
            f"[control] START  Tset = {t_set:.2f} C (low {low:g}, high {high:g}), "
            f"Kp = {kp:g} PWM/C, Ki = {ki:g} PWM/(C s), u_max = {self.run['u_max']:g}, "
            f"anti-windup {'on' if self.run['anti_windup'] else 'off'}, "
            f"logging to {self.csv_path}"
        )

    def stop_control(self, reason: str):
        if not self.control_on:
            return
        self.control_on = False
        self.send_line("SET PWM 0 DIR HEAT")

        if self.csv_file is not None:
            self.csv_file.close()
            self.csv_file = None
        print(f"[control] STOP ({reason}). q = {self.q:.1f} C s. Saved {self.csv_path}")

        # Back to manual mode at PWM 0. q keeps its value until "Zero integral".
        self.pwm_slider.blockSignals(True)
        self.pwm_slider.setValue(0)
        self.pwm_slider.blockSignals(False)
        self.pwm_edit.setText("0")
        self.dir_combo.blockSignals(True)
        self.dir_combo.setCurrentText("HEAT")
        self.dir_combo.blockSignals(False)
        for w in (self.low_spin, self.high_spin, self.kp_spin, self.ki_spin,
                  self.umax_spin, self.antiwindup_check, self.zero_button,
                  self.dir_combo, self.pwm_slider, self.pwm_edit):
            w.setEnabled(True)
        self.control_button.setText("Start control")
        for lbl, text in ((self.lbl_error, "e: ---"), (self.lbl_up, "u_P: ---"),
                          (self.lbl_ui, "u_I: ---"), (self.lbl_ureq, "u_req: ---"),
                          (self.lbl_u, "applied: ---"), (self.lbl_tset, "Tset: ---"),
                          (self.lbl_state, "")):
            lbl.setText(text)
        self.lbl_state.setStyleSheet("")

    def write_row(self, t, dt, T, t_set, step: PIStep, pwm_ard, hc, safety, control_on):
        r = self.run
        self.csv_file.write(
            f"{t:.2f},{dt:.3f},{T:.3f},{t_set:.2f},{r['kp']:g},{r['ki']:g},{r['u_max']:g},"
            f"{int(r['anti_windup'])},{step.e:.3f},{step.q:.3f},{step.u_p:.3f},"
            f"{step.u_i:.3f},{step.u_req:.3f},{step.u_applied:.3f},{int(step.held)},"
            f"{int(step.saturated)},{step.pwm},{step.direction},"
            f"{pwm_ard},{hc},{safety},{control_on}\n"
        )

    # ---------------- one line received ----------------
    @QtCore.Slot(str)
    def on_line(self, line: str):
        parsed = parse_measurement(line)
        if parsed is None:
            if "SAFETY" in line.upper():
                print(f"[arduino] {line}")
            return

        t = parsed["time_s"]
        T = parsed["temperature_C"]
        pwm_ard = parsed["pwm"]
        hc = parsed["heat_cool"]
        safety = parsed["safety_shutdown"]

        # dt = actual elapsed time since the previous accepted reading (Arduino
        # clock). A reset Arduino starts its clock again; then no integration, the
        # baseline restarts (a run file never mixes the two clocks), and the plots
        # continue to the right instead of drawing over the earlier readings.
        dt = 0.0 if self.t_prev is None else t - self.t_prev
        if dt < 0.0:
            print(f"[warning] Arduino time went back ({self.t_prev:.2f} -> {t:.2f} s): "
                  "the Arduino restarted; dt set to 0")
            dt = 0.0
            self.baseline.clear()
            if self.ts:
                self.t_offset = self.ts[-1] + 1.0 - t
        self.t_prev = t

        setpoint = error = u_p = u_i = u_req = float("nan")

        if self.control_on:
            r = self.run
            setpoint = self.active_setpoint()
            step = pi_control(setpoint, T, self.q, dt, r["kp"], r["ki"],
                              r["u_max"], r["anti_windup"])
            self.q = step.q
            error, u_p, u_i, u_req = step.e, step.u_p, step.u_i, step.u_req

            if safety:
                # The Arduino has already forced PWM to 0. Do not keep
                # commanding the TEC: leave control and require a restart.
                step = step._replace(u_applied=0.0, pwm=0, direction="HEAT")
            else:
                self.send_line(f"SET PWM {step.pwm} DIR {step.direction}", echo=False)

            self.write_row(t, dt, T, setpoint, step, pwm_ard, hc, safety, control_on=1)
            self.csv_file.flush()

            self.lbl_tset.setText(f"Tset: {setpoint:.2f} C")
            self.lbl_error.setText(f"e: {step.e:+.2f} C")
            self.lbl_up.setText(f"u_P: {step.u_p:+.1f}")
            self.lbl_ui.setText(f"u_I: {step.u_i:+.1f}")
            self.lbl_ureq.setText(f"u_req: {step.u_req:+.1f}")
            self.lbl_u.setText(f"applied: {step.u_applied:+.1f} (PWM {step.pwm} {step.direction})")
            flags = [name for name, on in (("SATURATED", step.saturated),
                                           ("INTEGRAL HELD", step.held)) if on]
            self.lbl_state.setText("  ".join(flags))
            self.lbl_state.setStyleSheet(
                "color: white; background-color: darkorange; font-weight: bold; padding: 2px;"
                if flags else "")
            print(
                f"t = {t:.2f} s, dt = {dt:.2f} s, T = {T:.2f} C, Tset = {setpoint:.2f} C, "
                f"e = {step.e:+.2f} C, q = {step.q:+.1f} C s, u_P = {step.u_p:+.1f}, "
                f"u_I = {step.u_i:+.1f}, u_req = {step.u_req:+.1f}, "
                f"sent PWM {step.pwm} {step.direction}"
                + ("  SATURATED" if step.saturated else "")
                + ("  HELD" if step.held else "")
                + f", Arduino PWM {pwm_ard} {'HEAT' if hc else 'COOL'}"
            )

            if safety:
                self.stop_control("Arduino safety shutdown")
        else:
            self.baseline.append((t, dt, T, pwm_ard, hc, safety))
            if safety and self.pwm_slider.value() != 0:
                # Keep the Module 4 behaviour: after a safety shutdown the
                # TEC stays off until the user sets a new PWM.
                self.pwm_slider.blockSignals(True)
                self.pwm_slider.setValue(0)
                self.pwm_slider.blockSignals(False)
                self.pwm_edit.setText("0")
            # Re-send the manual command once per reading so the Arduino's
            # command timeout does not stop an open-loop run.
            self.send_manual_command(echo=False)
            print(
                f"Temperature (C): {T:.2f}, Time (s): {t:.2f}, PWM: {pwm_ard}, "
                f"Heat/Cool: {hc}, Safety shutdown: {safety}"
            )

        # Store for the plots
        self.ts.append(t + self.t_offset)
        self.temps.append(T)
        self.setpoints.append(setpoint)
        self.errors.append(error)
        self.u_ps.append(u_p)
        self.u_is.append(u_i)
        self.u_reqs.append(u_req)
        self.pwm_signed.append(pwm_ard if hc == 1 else -pwm_ard)

        # Update the live labels
        self.lbl_temp.setText(f"T: {T:.2f} C")
        self.lbl_q.setText(f"q: {self.q:+.1f} C s")
        self.lbl_dt.setText(f"dt: {dt:.2f} s")
        self.lbl_pwm.setText(f"Arduino PWM: {pwm_ard}")
        self.lbl_dir.setText(f"Direction: {'HEAT' if hc == 1 else 'COOL'}")
        self.lbl_time.setText(f"Time: {t:.2f} s")
        if safety:
            self.lbl_safety.setText("SAFETY SHUTDOWN ACTIVE - PWM forced to 0")
            self.lbl_safety.setStyleSheet(
                "color: white; background-color: red; font-weight: bold; padding: 2px;"
            )
        else:
            self.lbl_safety.setText("Safety: OK")
            self.lbl_safety.setStyleSheet("")

    @QtCore.Slot(str)
    def on_serial_error(self, msg: str):
        print(f"[serial error] {msg}")
        if self.control_on:
            # Closes the run file; the zero command may not reach the Arduino.
            self.stop_control("serial error")

    # ---------------- plot updates ----------------
    def update_plots(self):
        # The whole session is drawn. pyqtgraph auto-ranges to all of it until
        # the user zooms or pans; "Show all" turns auto-range back on.
        if not self.ts:
            return
        xs = self.ts
        self.temp_curve.setData(xs, self.temps)
        self.setpoint_curve.setData(xs, self.setpoints)
        self.up_curve.setData(xs, self.u_ps)
        self.ui_curve.setData(xs, self.u_is)
        self.ureq_curve.setData(xs, self.u_reqs)
        self.pwm_curve.setData(xs, self.pwm_signed)
        self.error_curve.setData(xs, self.errors)

    def show_all(self):
        for p in self.plots:
            p.enableAutoRange()

    # ---------------- clean shutdown ----------------
    def closeEvent(self, event):
        self.stop_control("GUI closed")
        self.timer.stop()
        self.send_timer.stop()
        self.reader.stop()
        self.link.close()
        super().closeEvent(event)


# ------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------
def main():
    app = QtWidgets.QApplication(sys.argv)
    # Same white plot background as the Module 4 and 5 GUIs.
    pg.setConfigOptions(background="w", foreground="k")
    window = ControlGUI()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
