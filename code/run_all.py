"""Ejecuta toda la cadena: datos → E1–E3 → E4–E6 → E7 → comprobaciones SymPy.
Uso: python code/run_all.py [n_boot]   (n_boot por defecto 500)"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
n_boot = int(sys.argv[1]) if len(sys.argv) > 1 else 500
steps = [
    ("build_dataset", [sys.executable, "code/build_dataset.py"]),
    ("exp_timing (E1-E3)", [sys.executable, "code/exp_timing.py", str(n_boot)]),
    ("exp_quality (E4-E6)", [sys.executable, "code/exp_quality.py", str(n_boot)]),
    ("exp_model (E7)", [sys.executable, "code/exp_model.py", str(max(50, n_boot // 5))]),
    ("proofs_sympy", [sys.executable, "code/proofs_sympy.py"]),
]
log = (ROOT / "output" / "run.log").open("w")
for name, cmd in steps:
    t0 = time.time()
    print(f"== {name} ==", flush=True); log.write(f"\n== {name} ==\n"); log.flush()
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    out = "\n".join(l for l in (r.stdout + r.stderr).splitlines() if "Warning" not in l and "f = lambda" not in l)
    log.write(out + "\n"); log.flush()
    print(out[-3000:])
    print(f"   ({time.time() - t0:.0f} s, rc={r.returncode})", flush=True)
    if r.returncode != 0:
        print("ERROR en", name); sys.exit(1)
print("Cadena completa.")
