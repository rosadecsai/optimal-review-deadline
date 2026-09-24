"""E7: calibración del problema del editor y plazo óptimo T*.

Ingredientes empíricos (de E1–E5):
  * perfil de riesgo diario h(s) en días antes del plazo (invariante a T, ley de
    Parkinson), estimado con las cohortes 2024–2026 (T = 20–21) y extendido para
    s >= S_EXT con una ley de potencia h(s) = a (1+s)^-b;
  * riesgo diario constante mu tras el plazo (cola tardía);
  * penalización de calidad de los informes de última hora y tardíos (E4).
Con ellos se obtienen, para cualquier T, las cuotas pre/átomo/tardía, la demora
esperada y E q(T); el editor maximiza V(T). Validación fuera de muestra: el
perfil de 2024–2026 predice las cohortes 2020–2023 (T = 12–20).
También: comprobación de la forma cerrada (W de Lambert) y Monte Carlo del
modelo de parada óptima del revisor (invariancia a T y riesgo creciente).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
from model import EditorModel, ReviewerModel, calibrate_acceptance, tstar_closed_form, value_base  # noqa: E402
from utils import load_reviews, save_result, to_latex, cluster_bootstrap  # noqa: E402

CALIB_YEARS = [2024, 2025, 2026]
S_EXT = 15          # a partir de aquí se usa la ley de potencia (los recuentos empíricos son escasos)
ATOM = int(C.ATOM_DAYS)


# ----------------------------------------------------- perfil de riesgo -----
def hazard_profile(t_before: np.ndarray, s_max: int = 19) -> dict:
    s = np.floor(t_before).astype(int)
    h = {}
    for k in range(s_max, -1, -1):
        at_risk = (s <= k).sum()
        h[k] = float((s == k).sum() / at_risk) if at_risk else 0.0
    late = s[s < 0]
    # riesgo tardío constante: MLE geométrico sobre días de retraso (d = -s, d >= 1)
    d = -late
    mu = float(len(d) / d.sum()) if len(d) else 0.0
    ss = np.arange(2, S_EXT + 1)
    hh = np.array([max(h[k], 1e-6) for k in ss])
    b, loga = np.polyfit(np.log(1 + ss), np.log(hh), 1)
    return dict(h=h, mu=mu, pl_a=float(np.exp(loga)), pl_b=float(-b))


class DeliveryProfile:
    """Distribución de entrega para cualquier T a partir de un perfil de riesgo."""

    def __init__(self, prof: dict, max_late: int = 60):
        self.h = prof["h"]; self.mu = prof["mu"]; self.a = prof["pl_a"]; self.b = prof["pl_b"]
        self.max_late = max_late

    def hazard(self, s: int) -> float:
        if s < 0:
            return self.mu
        if s < S_EXT:
            return self.h.get(s, 0.0)
        return min(0.999, self.a * (1 + s) ** (-self.b))

    def delivery(self, T: int) -> dict:
        surv, out = 1.0, {}
        for s in range(T - 1, -self.max_late - 1, -1):
            hs = 1.0 if s == -self.max_late else self.hazard(s)
            out[s] = surv * hs
            surv *= 1 - hs
        return out

    def moments(self, T: int) -> dict:
        d = self.delivery(T)
        s = np.array(list(d)); p = np.array(list(d.values()))
        late = p[s < 0].sum()
        return dict(T=T, share_pre=float(p[s >= ATOM].sum()), share_atom=float(p[(s >= 0) & (s < ATOM)].sum()),
                    share_late=float(late), mean_s=float((s * p).sum()) + 0.5,
                    mean_lateness=float((-s[s < 0] * p[s < 0]).sum() / late - 0.5) if late > 0 else 0.0,
                    expected_delay=float(T - ((s * p).sum() + 0.5)))


def observed_moments(g: pd.DataFrame) -> dict:
    return dict(share_pre=float((g["phase"] == "pre").mean()), share_atom=float((g["phase"] == "atom").mean()),
                share_late=float(g["late"].mean()), mean_s=float(g["t_before"].mean()),
                mean_lateness=float((-g.loc[g["late"], "t_before"]).mean()),
                expected_delay=float((g["t_since_assign"]).mean()) if g["t_since_assign"].notna().any() else np.nan)


# ------------------------------------------------------------- escenarios ---
T_REF = 20.5   # ventana media de las cohortes de calibración


def scenarios(eta: float) -> list[dict]:
    """eta: elasticidad de compresión, Δ(T) = Δ·(T_ref/T)^eta (0 = penalización fija)."""
    abar, lam_a = calibrate_acceptance()
    S = []
    for delta in (0.0025, 0.005, 0.01):
        S.append(dict(name=f"Conference, δ={delta}, η={eta}", conference=True, delta=delta, L=0.0, abar=1.0, lam_a=1.0, eta=eta))
    for L in (0.1, 0.3, 0.6):
        for delta in (0.0025, 0.005, 0.01):
            S.append(dict(name=f"Journal, L={L}, δ={delta}, η={eta}", conference=False, delta=delta, L=L, abar=abar, lam_a=lam_a, eta=eta))
    return S


def value_curve(ed: EditorModel, Ts, moms, eta: float, Da: float, Dl: float, late_slope: float = 0.0):
    V = []
    for T, m in zip(Ts, moms):
        comp = (T_REF / T) ** eta
        ed.Delta_atom = Da * comp; ed.Delta_late = Dl * comp
        sl = m["share_late"]; delay = m["expected_delay"]
        if late_slope:
            sl2 = float(np.clip(sl + late_slope * (T - T_REF), 0.02, 0.8))
            delay += (sl2 - sl) * (m["mean_lateness"] + 5.5)   # masa trasladada de "pre" (≈5,5 días antes) a tardía
            sl = sl2
        V.append(ed.value(T, m["share_atom"], sl, delay))
    return np.array(V)


DEV_T = [7, 10, 14, 21, 28, 42]
ETA = 1.3          # elasticidad de compresión: penalización 2023 (T=12) frente a 2024-26 (T≈21) => (21/12)^η ≈ 2.1
LATE_SLOPE = -0.012  # pendiente de la cuota tardía sobre T entre cohortes (E3), variante de robustez


def e7(df: pd.DataFrame, Delta_atom: float, Delta_late: float, n_boot: int = 200):
    cal = df[df["year"].isin(CALIB_YEARS)]
    prof = hazard_profile(cal["t_before"].values)
    dp = DeliveryProfile(prof)
    print(f"perfil: mu={prof['mu']:.3f}; ley de potencia a={prof['pl_a']:.3f} b={prof['pl_b']:.3f}; h(0..5)=",
          [round(prof['h'][k], 3) for k in range(6)])

    # --- validación fuera de muestra: cohortes 2020-2023 (y 2024-26 dentro de muestra)
    val = []
    for y in sorted(C.YEARS):
        T = C.window_days(y, True)
        if T is None:
            continue
        g = df[df["year"] == y]
        obs = observed_moments(g); pred = dp.moments(int(T))
        val.append(dict(year=y, T=T, in_calibration=y in CALIB_YEARS,
                        **{f"obs_{k}": obs[k] for k in ["share_pre", "share_atom", "share_late", "mean_s", "mean_lateness"]},
                        **{f"pred_{k}": pred[k] for k in ["share_pre", "share_atom", "share_late", "mean_s", "mean_lateness"]}))
    val = pd.DataFrame(val)
    print(val.round(3).to_string())

    # --- curvas por T y escenarios
    Ts = list(range(3, 61))
    moms = [dp.moments(T) for T in Ts]
    curves = pd.DataFrame(moms)
    rows, Vcurves = [], {}
    for sc in scenarios(0.0) + scenarios(ETA):
        ed = EditorModel(delta=sc["delta"], L=sc["L"], abar=sc["abar"], lam_a=sc["lam_a"], Delta_atom=Delta_atom,
                         Delta_late=Delta_late, conference=sc["conference"])
        V = value_curve(ed, Ts, moms, sc["eta"], Delta_atom, Delta_late)
        i = int(np.argmax(V)); Tst = Ts[i]
        # T* con IC bootstrap sobre el perfil (remuestreo de manuscritos de las cohortes de calibración)
        def stat(dd):
            p2 = DeliveryProfile(hazard_profile(dd["t_before"].values))
            mm = [p2.moments(T) for T in Ts]
            VV = value_curve(ed, Ts, mm, sc["eta"], Delta_atom, Delta_late)
            return [Ts[int(np.argmax(VV))]]
        _, lo, hi = cluster_bootstrap(cal, stat, n_boot=n_boot, seed=C.SEED + 7)
        Vls = value_curve(ed, Ts, moms, sc["eta"], Delta_atom, Delta_late, late_slope=LATE_SLOPE)
        r = dict(scenario=sc["name"], conference=sc["conference"], delta=sc["delta"], L=sc["L"], eta=sc["eta"], T_star=Tst, T_star_lo=float(lo[0]), T_star_hi=float(hi[0]),
                 T_star_late_slope=Ts[int(np.argmax(Vls))],
                 V_star=float(V[i]), a_Tstar=float(ed.a(Tst)), share_atom_Tstar=moms[i]["share_atom"], share_late_Tstar=moms[i]["share_late"],
                 expected_delay_Tstar=moms[i]["expected_delay"])
        for Td in DEV_T:
            r[f"loss_T{Td}"] = float(100 * (V[i] - V[Ts.index(Td)]) / ed.qbar)  # % de q̄
        rows.append(r); Vcurves[sc["name"]] = V
        print(f"{sc['name']:34s} T*={Tst:3d} [{lo[0]:.0f},{hi[0]:.0f}] (late-slope variant {r['T_star_late_slope']})  V*={V[i]:.3f}  loss@7={r['loss_T7']:.2f}% @14={r['loss_T14']:.2f}% @28={r['loss_T28']:.2f}% @42={r['loss_T42']:.2f}%")
    tab6 = pd.DataFrame(rows)

    # --- sensibilidad a la penalización de calidad Δ (escenario central de revista)
    abar, lam_a = calibrate_acceptance()
    sens = []
    for Da in (0.05, 0.085, 0.15, 0.3):
        for delta in (0.0025, 0.005, 0.01):
            for la in (0.05, lam_a, 0.2):
                for eta in (0.0, ETA):
                    ed = EditorModel(delta=delta, L=0.3, abar=abar, lam_a=la, conference=False)
                    V = value_curve(ed, Ts, moms, eta, Da, Da)
                    edc = EditorModel(delta=delta, conference=True)
                    Vc = value_curve(edc, Ts, moms, eta, Da, Da)
                    sens.append(dict(Delta=Da, delta=delta, lam_a=round(la, 4), eta=eta, T_star_journal=Ts[int(np.argmax(V))],
                                     T_star_conference=Ts[int(np.argmax(Vc))]))
    sens = pd.DataFrame(sens)

    # --- figura 4
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    ax = axes[0]
    ax.plot(curves["T"], 100 * curves["share_atom"], label="last 48 h (atom)")
    ax.plot(curves["T"], 100 * curves["share_late"], label="late")
    ax.plot(curves["T"], 100 * curves["share_pre"], label="earlier")
    for _, r in val.iterrows():
        ax.plot(r["T"], 100 * r["obs_share_atom"], "o", color="C0", ms=4); ax.plot(r["T"], 100 * r["obs_share_late"], "s", color="C1", ms=4)
    ax.set_xlabel("Review window T (days)"); ax.set_ylabel("Share of reviews (%)"); ax.set_title("(a) Predicted shares vs T (dots: observed cohorts)")
    ax.legend(fontsize=8, frameon=False)
    ax = axes[1]
    ax.plot(curves["T"], curves["expected_delay"], label="expected delivery (days after assignment)")
    ax.plot(curves["T"], curves["T"], ls="--", color="grey", lw=0.8, label="T")
    ax.set_xlabel("Review window T (days)"); ax.set_ylabel("Days"); ax.set_title("(b) Expected delivery time vs T")
    ax.legend(fontsize=8, frameon=False)
    ax = axes[2]
    for name, V in Vcurves.items():
        if name in (f"Conference, δ=0.005, η=0.0", f"Conference, δ=0.005, η={ETA}", "Journal, L=0.3, δ=0.005, η=0.0", f"Journal, L=0.3, δ=0.005, η={ETA}", "Journal, L=0.3, δ=0.01, η=0.0"):
            ax.plot(Ts, V, label=name)
    ax.set_xlabel("Review window T (days)"); ax.set_ylabel("Editor's expected value V(T)"); ax.set_title("(c) V(T) and T* by scenario")
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout(); fig.savefig(C.FIG_DIR / "fig4_optimal_deadline.pdf"); fig.savefig(C.FIG_DIR / "fig4_optimal_deadline.png", dpi=200); plt.close(fig)

    # --- forma cerrada (caso base) frente al óptimo numérico
    cf = []
    for kw in [dict(), dict(delta=0.01), dict(L=0.6), dict(lam_a=0.2, m=2.0), dict(qbar=1.5, delta=0.0025, m=1.5)]:
        Tg = np.linspace(0.5, 200, 200000)
        Tn = float(Tg[np.argmax(value_base(Tg, **kw))]); Tc = tstar_closed_form(**kw)
        cf.append(dict(params=str(kw), T_closed_form=Tc, T_numeric=Tn, abs_diff=abs(Tc - Tn)))
    cf = pd.DataFrame(cf); print(cf.round(4).to_string())

    return prof, val, curves, tab6, sens, cf


# ------------------------------------------ Monte Carlo del modelo del revisor
def reviewer_monte_carlo(n: int = 20000):
    """Tipos heterogéneos (sigma, beta); misma población para T = 12, 16, 21, 28, 45.
    Muestra que la distribución en días antes del plazo es invariante a T salvo
    truncamiento y que el riesgo crece hacia el plazo."""
    g = np.random.Generator(np.random.PCG64(C.SEED))
    types = [(0.4, 0.5, 0.35), (0.8, 0.7, 0.35), (1.2, 0.9, 0.2), (0.6, 1.0, 0.1)]
    out = {}
    for T in (12, 16, 21, 28, 45):
        dens = {}
        for sig, beta, wgt in types:
            m = ReviewerModel(sigma=sig, beta=beta, late_pen=0.25, late_growth=0.02)
            d = m.delivery(T)
            for s, p in d.items():
                dens[s] = dens.get(s, 0.0) + wgt * p
        out[T] = dens
    rows = []
    for T, dens in out.items():
        s = np.array(list(dens)); p = np.array(list(dens.values()))
        rows.append(dict(T=T, share_atom=float(p[(s >= 0) & (s < ATOM)].sum()), share_late=float(p[s < 0].sum()),
                         mean_s=float((s * p).sum()), q50_s=float(s[np.searchsorted(np.cumsum(p[np.argsort(-s)]), 0.5)] if True else np.nan)))
    fig, ax = plt.subplots(figsize=(6, 4))
    for T, dens in out.items():
        s = np.array(sorted(dens)); p = np.array([dens[k] for k in s])
        ax.plot(s, 100 * p, label=f"T={T}")
    ax.invert_xaxis(); ax.set_xlim(30, -10); ax.set_xlabel("Days before the deadline"); ax.set_ylabel("Share of reviews (%)")
    ax.set_title("Optimal-stopping model: delivery distribution for different T"); ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(C.FIG_DIR / "figS1_model_invariance.pdf"); fig.savefig(C.FIG_DIR / "figS1_model_invariance.png", dpi=200); plt.close(fig)
    return pd.DataFrame(rows)


def main(n_boot: int = 200):
    df = load_reviews()
    import json
    res = json.loads((C.ROOT / "output" / "results.json").read_text())
    e4 = pd.DataFrame(res["E4"]["pooled"])
    Da = float(-e4[(e4["outcome"] == "log_words") & (e4["bin"] == "last24h")]["coef"].iloc[0])   # ≈ 0.085
    Dl = float(-e4[(e4["outcome"] == "log_words") & (e4["bin"] == "late")]["coef"].iloc[0])      # ≈ 0.095
    print(f"Δ_atom={Da:.3f} Δ_late={Dl:.3f} (penalización de longitud, efectos fijos de manuscrito)")
    prof, val, curves, tab6, sens, cf = e7(df, Da, Dl, n_boot)
    to_latex(val, C.TAB_DIR / "table_e7_validation", "Out-of-sample prediction of cohort timing from the 2024–2026 hazard profile", "tab:validation", "%.3f")
    to_latex(curves, C.TAB_DIR / "table_e7_curves", "Predicted timing moments by window", "tab:curves", "%.3f")
    to_latex(tab6, C.TAB_DIR / "table6_tstar", "Optimal review window by scenario", "tab:tstar", "%.2f")
    to_latex(sens, C.TAB_DIR / "table6_sensitivity", "Sensitivity of T* to the quality penalty", "tab:sens", "%.3f")
    print(sens.pivot_table(index=["Delta", "eta"], columns=["lam_a", "delta"], values="T_star_journal").to_string())
    print(sens.pivot_table(index=["Delta", "eta"], columns="delta", values="T_star_conference").to_string())
    mc = reviewer_monte_carlo()
    mc.to_csv(C.TAB_DIR / "tableS_model_mc.csv", index=False); print(mc.round(3).to_string())
    save_result("E7", dict(profile=dict(mu=prof["mu"], pl_a=prof["pl_a"], pl_b=prof["pl_b"], h={str(k): v for k, v in prof["h"].items()}),
                           Delta_atom=Da, Delta_late=Dl, validation=val.to_dict(orient="records"), tstar=tab6.to_dict(orient="records"),
                           sensitivity=sens.to_dict(orient="records"), closed_form=cf.to_dict(orient="records"), model_mc=mc.to_dict(orient="records")))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 200)
