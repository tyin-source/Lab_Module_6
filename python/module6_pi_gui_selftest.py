# Phys 39 Module 6, Part 7: self-test of the P/PI GUI without hardware
#
# Runs module6_pi_control_gui.py off screen against a simulated Arduino: a one-lump
# block (tau = 70 s, chi_h = 0.485, chi_c = 0.172 C/PWM) that prints one line per
# second in the Arduino format and applies the last command it received. Checks:
#   1. the Part 6 windup case (Kp 8, Ki 0.36, u_max 30, 19.84 -> 30 C) with
#      anti-windup on and off; the Part 6 model gives T_max = 30.20 and 33.27 C
#   2. that an Arduino safety-shutdown line stops control and sends PWM 0
# Run files go to a temporary folder, not data/module_06/.
#
# Usage (from the repository root):
#   python python/module6_pi_gui_selftest.py

import contextlib
import csv
import io
import os
import sys
import tempfile
import time
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6 import QtWidgets
import pyqtgraph as pg

import module6_pi_control_gui as gui_module

CHI_H, CHI_C, TAU = 0.485, 0.172, 70.0
T_AMB = 19.84


class SimulatedArduino:
    """Fake serial link and one-lump block in one object."""

    def __init__(self, t_amb=T_AMB):
        self.pwm, self.heat, self.sent = 0, True, []
        self.T, self.t, self.t_amb = t_amb, 0.0, t_amb

    # serial link interface used by the GUI
    def readline(self):
        time.sleep(0.05)
        return b""

    def write_line(self, text):
        self.sent.append(text)
        words = text.split()
        self.pwm, self.heat = int(words[2]), words[4] == "HEAT"

    def close(self):
        pass

    # one ~1 s average with the last command applied, then one printed line
    def line(self, safety=0):
        u = self.pwm if self.heat else -self.pwm
        chi = CHI_H if u >= 0 else CHI_C
        for _ in range(100):
            self.T += 0.01 / TAU * (self.t_amb + chi * u - self.T)
        self.t += 1.0
        return (f"Temperature (C): {self.T:.2f}, Time (s): {self.t:.2f}, PWM: {self.pwm}, "
                f"Direction input: {int(self.heat)}, Active PWM pin: {9 if self.heat else 10}, "
                f"Heat/Cool: {int(self.heat)}, Safety shutdown: {safety}")


def feed(app, gui, arduino, n):
    with contextlib.redirect_stdout(io.StringIO()):
        for _ in range(n):
            gui.on_line(arduino.line())
            app.processEvents()


def windup_run(app, anti_windup):
    arduino = SimulatedArduino()
    gui = gui_module.ControlGUI(link=arduino)
    feed(app, gui, arduino, 35)                      # baseline at PWM 0
    gui.low_spin.setValue(19.0)
    gui.high_spin.setValue(30.0)
    gui.kp_spin.setValue(8.0)
    gui.ki_spin.setValue(0.36)
    gui.umax_spin.setValue(30)
    gui.antiwindup_check.setChecked(anti_windup)
    gui.high_button.click()
    with contextlib.redirect_stdout(io.StringIO()):
        gui.control_button.click()
    feed(app, gui, arduino, 400)
    path = gui.csv_path
    with contextlib.redirect_stdout(io.StringIO()):
        gui.close()
    rows = [r for r in csv.DictReader(open(path)) if r["control_on"] == "1"]
    T = [float(r["temperature_C"]) for r in rows]
    t0 = float(rows[0]["time_s"])
    k = max(range(len(T)), key=T.__getitem__)
    print(f"  anti-windup {'on ' if anti_windup else 'off'}: T_max {T[k]:.2f} C at "
          f"{float(rows[k]['time_s']) - t0:.0f} s, largest u_I "
          f"{max(float(r['u_I']) for r in rows):.1f} PWM, steps with q held "
          f"{sum(r['integral_held'] == '1' for r in rows)}, largest pwm_cmd "
          f"{max(int(r['pwm_cmd']) for r in rows)}   ({path.name})")
    return T[k]


def safety_run(app):
    arduino = SimulatedArduino()
    gui = gui_module.ControlGUI(link=arduino)
    feed(app, gui, arduino, 5)
    gui.kp_spin.setValue(4.0)
    gui.ki_spin.setValue(0.1)
    with contextlib.redirect_stdout(io.StringIO()):
        gui.control_button.click()
    feed(app, gui, arduino, 5)
    n = len(arduino.sent)
    with contextlib.redirect_stdout(io.StringIO()):
        gui.on_line(arduino.line(safety=1))
    ok = not gui.control_on and arduino.sent[n:] == ["SET PWM 0 DIR HEAT"]
    print(f"  safety-shutdown line: control stopped and PWM 0 sent: {'yes' if ok else 'NO'}")
    with contextlib.redirect_stdout(io.StringIO()):
        gui.close()
    return ok


def main():
    app = QtWidgets.QApplication(sys.argv)
    pg.setConfigOptions(background="w", foreground="k")
    with tempfile.TemporaryDirectory() as tmp:
        gui_module.DATA_DIR = Path(tmp)
        print("Windup case of Part 6 (model: T_max 30.20 C with anti-windup, 33.27 C without)")
        on = windup_run(app, True)
        off = windup_run(app, False)
        ok = safety_run(app)
    passed = abs(on - 30.20) < 0.1 and abs(off - 33.27) < 0.1 and ok
    print("PASS" if passed else "FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
