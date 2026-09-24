"""Comprobaciones simbólicas (SymPy) y numéricas de las proposiciones del modelo.

P1  Invariancia (parada óptima con plazo): los umbrales dependen solo de los días
    que faltan, luego el riesgo h(s) es el mismo para cualquier T (comprobación
    numérica con el modelo del revisor).
P2  Techo de calidad: E[q|T] = q̄ − Δ[(1−π)e^{−λT} + πφ] es creciente, cóncava y
    tiende a q̄ − Δπφ < q̄.
P3  Existencia/unicidad y forma cerrada: en el caso base V(T) = a(T)(q̄ − δ(T−m))
    − (1−a(T))L con a(T) = ā(1−e^{−λ_a T}), la CPO tiene solución única
    T* = [A/δ − W(e^{A/δ})]/λ_a con A/δ = 1 + λ_a(q̄ + L + δm)/δ, y V''(T*) < 0.
P4  Estática comparativa: dT*/dδ < 0, dT*/dL > 0, dT*/dq̄ > 0, dT*/dλ_a < 0 (signo
    por diferenciación implícita, comprobado en una rejilla amplia).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from model import ReviewerModel, tstar_closed_form  # noqa: E402

checks = []


def check(name, ok):
    checks.append((name, bool(ok)))
    print(("OK  " if ok else "FAIL") + " " + name)


# --- P1: invariancia numérica del riesgo ---------------------------------------
m = ReviewerModel(sigma=0.8, beta=0.7, late_pen=0.25, late_growth=0.02)
def haz_from_delivery(T):
    d = m.delivery(T); surv = 1.0; h = {}
    for s in sorted(d, reverse=True):
        h[s] = d[s] / surv if surv > 1e-12 else np.nan
        surv -= d[s]
    return h
h12, h45 = haz_from_delivery(12), haz_from_delivery(45)
common = [s for s in h12 if s in h45 and s < 12]
check("P1 riesgo h(s) idéntico para T=12 y T=45 en los días comunes",
      max(abs(h12[s] - h45[s]) for s in common) < 1e-9)
check("P1 riesgo creciente hacia el plazo (s de 11 a 0)", all(h12[s] >= h12[s + 1] for s in range(0, 11)))

# --- P2: techo de calidad -------------------------------------------------------
T, q, D, pi, lam, phi = sp.symbols("T qbar Delta pi lambda phi", positive=True)
Eq = q - D * ((1 - pi) * sp.exp(-lam * T) + pi * phi)
dEq = sp.simplify(sp.diff(Eq, T)); d2Eq = sp.simplify(sp.diff(Eq, T, 2))
check("P2 dE[q]/dT = Δ(1−π)λe^{−λT} > 0", sp.simplify(dEq - D * (1 - pi) * lam * sp.exp(-lam * T)) == 0)
check("P2 d²E[q]/dT² = −Δ(1−π)λ²e^{−λT} < 0", sp.simplify(d2Eq + D * (1 - pi) * lam ** 2 * sp.exp(-lam * T)) == 0)
check("P2 límite T→∞: q̄ − Δπφ", sp.simplify(sp.limit(Eq, T, sp.oo) - (q - D * pi * phi)) == 0)

# --- P3: forma cerrada ----------------------------------------------------------
L, delta, la, abar, mm = sp.symbols("L delta lambda_a abar m", positive=True)
a = abar * (1 - sp.exp(-la * T))
V = a * (q - delta * (T - mm)) - (1 - a) * L
dV = sp.diff(V, T)
# CPO: e^{-λT}[λ(q̄+L+δm) + δ − λδT] = δ   (tras dividir por ā)
foc = sp.simplify(dV / abar - (sp.exp(-la * T) * (la * (q + L + delta * mm) + delta - la * delta * T) - delta))
check("P3 la CPO se reduce a e^{−λT}[λ(q̄+L+δm)+δ−λδT] = δ", foc == 0)
Ad = 1 + la * (q + L + delta * mm) / delta
Tstar = (Ad - sp.LambertW(sp.exp(Ad))) / la
# sustituir: con u = A/δ − λT* = W(e^{A/δ}) se cumple u e^{u} = e^{A/δ}, luego e^{−λT*}(λδ/δ)(...)...
expr = sp.exp(-la * Tstar) * (la * (q + L + delta * mm) + delta - la * delta * Tstar) - delta
# comprobación numérica en varios puntos (la simplificación simbólica de W es frágil)
pts = [dict(q=1, L=0.3, delta=0.005, la=0.1, mm=1), dict(q=1.5, L=0.6, delta=0.01, la=0.2, mm=2), dict(q=1, L=0.1, delta=0.0025, la=0.05, mm=0.5)]
ok = True
for p in pts:
    val = float(expr.subs({q: p["q"], L: p["L"], delta: p["delta"], la: p["la"], mm: p["mm"]}))
    ok &= abs(val) < 1e-9
check("P3 T* de forma cerrada satisface la CPO (3 puntos numéricos, |resid|<1e-9)", ok)
d2V = sp.diff(V, T, 2)
ok = True
for p in pts:
    sub = {q: p["q"], L: p["L"], delta: p["delta"], la: p["la"], mm: p["mm"], abar: 0.8}
    Ts = float(Tstar.subs(sub)); ok &= float(d2V.subs(sub).subs(T, Ts)) < 0
    # unicidad: dV cambia de signo una sola vez en (0, 500)
    grid = np.linspace(0.01, 500, 20000); vals = np.array([float(dV.subs(sub).subs(T, t)) for t in grid[::200]])
    ok &= int((np.diff(np.sign(vals)) != 0).sum()) == 1
check("P3 V''(T*) < 0 y dV cambia de signo una sola vez (unicidad)", ok)
# equivalencia con la implementación numérica en model.py
ok = all(abs(float(Tstar.subs({q: p["q"], L: p["L"], delta: p["delta"], la: p["la"], mm: p["mm"]})) -
             tstar_closed_form(qbar=p["q"], L=p["L"], delta=p["delta"], lam_a=p["la"], m=p["mm"])) < 1e-8 for p in pts)
check("P3 forma cerrada de model.py coincide con la simbólica", ok)

# --- P4: estática comparativa por diferenciación implícita ---------------------
G = sp.exp(-la * T) * (la * (q + L + delta * mm) + delta - la * delta * T) - delta   # G(T*, ·) = 0
GT = sp.diff(G, T)
sgn = {}
rng = np.random.default_rng(0)
for name, par in [("delta", delta), ("L", L), ("qbar", q), ("lambda_a", la)]:
    dTdp = -sp.diff(G, par) / GT
    signs = []
    for _ in range(300):
        p = dict(q=rng.uniform(0.5, 3), L=rng.uniform(0.0, 2), delta=rng.uniform(0.001, 0.05), la=rng.uniform(0.02, 0.5), mm=rng.uniform(0, 5))
        sub = {q: p["q"], L: p["L"], delta: p["delta"], la: p["la"], mm: p["mm"]}
        Ts = float(Tstar.subs(sub))
        if Ts <= 0:
            continue
        signs.append(np.sign(float(dTdp.subs(sub).subs(T, Ts))))
    sgn[name] = (min(signs), max(signs), len(signs))
check("P4 dT*/dδ < 0 en 300 puntos aleatorios", sgn["delta"][1] < 0)
check("P4 dT*/dL > 0", sgn["L"][0] > 0)
check("P4 dT*/dq̄ > 0", sgn["qbar"][0] > 0)
check("P4 dT*/dλ_a < 0", sgn["lambda_a"][1] < 0)
print("signos (min, max, n):", sgn)

n_ok = sum(ok for _, ok in checks)
print(f"\n{n_ok}/{len(checks)} comprobaciones superadas")
