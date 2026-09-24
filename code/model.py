"""Modelo teórico del plazo óptimo de revisión.

(A) Revisor: parada óptima con plazo (coste diario aleatorio, sesgo presente
    cuasi-hiperbólico, penalización por retraso). Genera la distribución de
    entrega en días antes del plazo, invariante a T salvo truncamiento (ley de
    Parkinson), con átomo en el plazo y cola tardía.
(B) Editor: V(T) = a(T)·[E q(T) − δ·E τ(T)] − (1 − a(T))·L, con a(T) creciente y
    cóncava (aceptación de la invitación), E q(T) creciente y cóncava con techo,
    E τ(T) ≈ T − m (anclaje al plazo). Forma cerrada de T* con la función W de
    Lambert en el caso base; óptimo numérico en el caso general.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import lambertw
from scipy.stats import lognorm


# ------------------------------------------------------------- (A) revisor ---
@dataclass
class ReviewerModel:
    """Coste diario c ~ LogNormal(mu, sigma) (media normalizada a 1 si mu=-sigma^2/2).
    Cada día el revisor hace el informe si su coste no supera el umbral (valor de
    esperar); con sesgo presente beta<1 el umbral es beta veces la continuación.
    Tras el plazo cada día de retraso cuesta `late_pen` (creciente) y hay un
    horizonte máximo de retraso `max_late`; en el último día se entrega seguro."""
    sigma: float = 1.0
    beta: float = 1.0
    late_pen: float = 0.5        # penalización por día de retraso (en unidades de coste medio)
    late_growth: float = 0.0     # crecimiento de la penalización con los días de retraso
    max_late: int = 30
    sophisticated: bool = True
    n_grid: int = 4000
    _cache: dict = field(default_factory=dict, repr=False)

    def cost_dist(self):
        return lognorm(s=self.sigma, scale=np.exp(-self.sigma ** 2 / 2))

    def thresholds(self):
        """Umbrales de aceptación w[k] para k = días hasta el final del horizonte
        (k=0: último día, entrega segura). Devuelve (w, W) con W la continuación."""
        D = self.cost_dist()
        # rejilla de cuantiles para esperanzas truncadas
        u = (np.arange(self.n_grid) + 0.5) / self.n_grid
        c = D.ppf(u)
        H = self.max_late
        # índice k cuenta días hasta el último día admisible (deadline + H)
        # penalización por retraso: día con s días antes del plazo (s<0 => retraso -s)
        def pen(s):
            d = max(0, -s)
            return self.late_pen * d * (1 + self.late_growth * d) if d > 0 else 0.0
        # calcular para un T grande y guardar por s (invariante a T): s = T..-H
        S_MAX = 200
        w = {}   # umbral efectivo (comparación con coste de hoy) por s
        W = {}   # coste esperado de continuar a partir de mañana (valor de espera), evaluado hoy
        Wtrue = {}  # continuación con el comportamiento real (para el sofisticado)
        # último día: s = -H → entrega segura
        s = -H
        Wtrue[s] = float(np.mean(c)) + pen(s)
        w[s] = np.inf
        for s in range(-H + 1, S_MAX + 1):
            nxt = Wtrue[s - 1]          # coste esperado si se pospone (comportamiento real mañana)
            # umbral hoy: hacer si c + pen(s) <= beta*nxt  (sesgo presente: el futuro pesa beta)
            thr = self.beta * nxt - pen(s)
            w[s] = thr
            # coste esperado hoy visto desde ayer (comportamiento real): E[(c+pen) 1{c<=thr}] + P(c>thr)*nxt
            done = c <= thr
            if self.sophisticated or self.beta == 1.0:
                Wtrue[s] = float(np.mean(np.where(done, c + pen(s), nxt)))
            else:
                # ingenuo: cree que mañana actuará sin sesgo (umbral con beta=1) → continuación percibida
                # se calcula aparte; aquí Wtrue sigue el comportamiento real (umbral con beta)
                Wtrue[s] = float(np.mean(np.where(done, c + pen(s), nxt)))
        if not self.sophisticated and self.beta < 1.0:
            # ingenuo: umbral usa la continuación de un yo futuro sin sesgo
            W1 = {}
            W1[-H] = float(np.mean(c)) + pen(-H)
            for s in range(-H + 1, S_MAX + 1):
                thr1 = W1[s - 1] - pen(s)
                W1[s] = float(np.mean(np.where(c <= thr1, c + pen(s), W1[s - 1])))
            for s in range(-H + 1, S_MAX + 1):
                w[s] = self.beta * W1[s - 1] - pen(s)
        return w, Wtrue

    def hazard(self):
        """h(s) = P(entrega en el día s | no entregado antes), s = días antes del plazo."""
        if "h" in self._cache:
            return self._cache["h"]
        w, _ = self.thresholds()
        D = self.cost_dist()
        h = {s: (1.0 if np.isinf(v) else float(D.cdf(v)) if v > 0 else 0.0) for s, v in w.items()}
        self._cache["h"] = h
        return h

    def delivery(self, T: int) -> dict:
        """Distribución de entrega para ventana T: dict s -> prob., s = T-1..-max_late.
        El día de asignación es s = T-1 (T días disponibles incluido el del plazo, s=0)."""
        h = self.hazard()
        surv, out = 1.0, {}
        for s in range(T - 1, -self.max_late - 1, -1):
            p = surv * h[s]
            out[s] = p
            surv *= 1 - h[s]
        return out

    def moments(self, T: int, atom_days: int = 2) -> dict:
        d = self.delivery(T)
        s = np.array(list(d.keys())); p = np.array(list(d.values()))
        pre = p[s >= atom_days].sum(); atom = p[(s >= 0) & (s < atom_days)].sum(); late = p[s < 0].sum()
        mean_s = float((s * p).sum())
        mean_late = float((-s[s < 0] * p[s < 0]).sum() / late) if late > 0 else 0.0
        # calidad: coste esperado del día de entrega (mayor coste = menos horas disponibles)
        return dict(T=T, share_pre=pre, share_atom=atom, share_late=late, mean_s=mean_s, mean_lateness=mean_late,
                    expected_delay=T - 1 - mean_s)


def empirical_density(t_before: np.ndarray, s_max: int = 25, s_min: int = -12) -> tuple[np.ndarray, np.ndarray]:
    s = np.floor(t_before).astype(int)
    ks = np.arange(s_min, s_max + 1)
    dens = np.array([(s == k).mean() for k in ks])
    return ks, dens


def fit_reviewer_model(t_before: np.ndarray, T: int, sophisticated: bool = True, s_max: int = 25, s_min: int = -12,
                       fit_beta: bool = True) -> tuple[ReviewerModel, dict]:
    """Ajusta (sigma, beta, late_pen, late_growth) a la densidad empírica diaria de días
    antes del plazo (mínimos cuadrados sobre raíz de la densidad) para una ventana T."""
    ks, dens = empirical_density(t_before, s_max, s_min)

    def loss(theta):
        sigma, beta, lp, lg = theta
        if sigma <= 0.05 or beta <= 0.05 or beta > 1.0 or lp <= 0 or lg < 0:
            return 1e6
        m = ReviewerModel(sigma=sigma, beta=beta, late_pen=lp, late_growth=lg, sophisticated=sophisticated)
        d = m.delivery(T)
        pred = np.array([d.get(k, 0.0) for k in ks])
        return float(np.sum((np.sqrt(pred) - np.sqrt(dens)) ** 2))

    best = None
    for x0 in ([1.0, 0.8, 0.3, 0.05], [1.5, 0.6, 0.2, 0.02], [0.7, 0.9, 0.5, 0.1], [2.0, 0.7, 0.1, 0.01]):
        if not fit_beta:
            x0[1] = 1.0
        r = minimize(loss, x0, method="Nelder-Mead", options=dict(maxiter=600, xatol=1e-3, fatol=1e-7))
        if best is None or r.fun < best.fun:
            best = r
    sigma, beta, lp, lg = best.x
    m = ReviewerModel(sigma=sigma, beta=beta, late_pen=lp, late_growth=lg, sophisticated=sophisticated)
    d = m.delivery(T)
    pred = np.array([d.get(k, 0.0) for k in ks])
    return m, dict(loss=best.fun, ks=ks, dens=dens, pred=pred, params=dict(sigma=sigma, beta=beta, late_pen=lp, late_growth=lg))


# --------------------------------------------------------------- (B) editor ---
@dataclass
class EditorModel:
    qbar: float = 1.0        # valor de un informe completo
    delta: float = 0.005     # coste de demora por día (fracción de qbar)
    L: float = 0.3           # pérdida si el revisor rechaza la invitación
    abar: float = 0.85       # aceptación máxima
    lam_a: float = 0.10      # velocidad con que la aceptación crece con T
    Delta_atom: float = 0.12 # pérdida de calidad de un informe de última hora (fracción de qbar)
    Delta_late: float = 0.12 # pérdida de calidad de un informe tardío
    conference: bool = False # a(T)=1, L=0 (revisores precomprometidos)

    def a(self, T):
        if self.conference:
            return np.ones_like(np.asarray(T, dtype=float))
        return self.abar * (1 - np.exp(-self.lam_a * np.asarray(T, dtype=float)))

    def value(self, T, share_atom, share_late, expected_delay):
        Eq = self.qbar - self.Delta_atom * share_atom - self.Delta_late * share_late
        a = self.a(T)
        L = 0.0 if self.conference else self.L
        return a * (Eq - self.delta * expected_delay) - (1 - a) * L

    def curve(self, reviewer: ReviewerModel, Ts=range(3, 61)) -> dict:
        Ts = list(Ts)
        mom = [reviewer.moments(T) for T in Ts]
        V = np.array([self.value(T, m["share_atom"], m["share_late"], m["expected_delay"]) for T, m in zip(Ts, mom)])
        i = int(np.argmax(V))
        return dict(T=np.array(Ts), V=V, T_star=Ts[i], V_star=V[i], moments=mom)


def tstar_closed_form(qbar=1.0, L=0.3, delta=0.005, lam_a=0.10, m=1.0):
    """Caso base sin compresión y entrega en T − m: FOC resuelta con W de Lambert:
       A/δ = 1 + λ_a (q̄ + L + δ m)/δ ;  T* = [A/δ − W(e^{A/δ})]/λ_a."""
    Ad = 1 + lam_a * (qbar + L + delta * m) / delta
    return float((Ad - lambertw(np.exp(Ad)).real) / lam_a)


def value_base(T, qbar=1.0, L=0.3, delta=0.005, lam_a=0.10, abar=0.85, m=1.0):
    a = abar * (1 - np.exp(-lam_a * T))
    return a * (qbar - delta * (T - m)) - (1 - a) * L


def calibrate_acceptance(a28=0.641, a45=0.676):
    """a(T) = ā(1 − e^{−λ_a T}) que pasa por los dos plazos del RCT de Chetty et al. (2014)."""
    from scipy.optimize import brentq
    f = lambda lam: (1 - np.exp(-45 * lam)) / (1 - np.exp(-28 * lam)) - a45 / a28
    lam = brentq(f, 1e-4, 1.0)
    abar = a28 / (1 - np.exp(-28 * lam))
    return abar, lam
