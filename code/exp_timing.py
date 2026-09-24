"""E1 distribución de entrega respecto al plazo; E2 mezcla (pre-plazo anclado al
plazo + átomo del día del plazo + cola tardía) frente al modelo de riesgo
constante desde la asignación; E3 ley de Parkinson entre ventanas anuales."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
from utils import cluster_bootstrap, load_reviews, save_result, to_latex, rng  # noqa: E402

S0 = C.ATOM_DAYS  # inicio del átomo del plazo (últimas 48 h) en días antes de las 23:59 AoE
YEARS = sorted(C.YEARS)


# ---------------------------------------------------------------- E1 ---------
def e1_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for y, g in df.groupby("year"):
        late = g[g["late"]]
        lateness = -late["t_before"]
        daily = np.floor(-g["t_before"]).astype(int).value_counts().sort_index()  # días desde el plazo
        peak = daily.max()
        # picos secundarios tras el plazo: máximos locales > 25 % del pico, a >= 2 días
        post = daily[daily.index >= 2]
        sec = [int(k) for k in post.index if post[k] > 0.25 * peak and post[k] >= post.get(k - 1, 0) and post[k] >= post.get(k + 1, 0)]
        rows.append(dict(
            year=y, window_official=C.window_days(y), window_effective=C.window_days(y, True),
            n_papers=g["forum"].nunique(), n_reviews=len(g),
            median_days_before=float(g["t_before"].median()),
            share_pre=float((g["phase"] == "pre").mean()), share_deadline_date=float((g["phase"] == "atom").mean()),
            share_last24h=float((g["bin"] == "last24h").mean()), share_late=float(g["late"].mean()),
            mean_lateness_days=float(lateness.mean()), median_lateness_days=float(lateness.median()),
            p90_lateness_days=float(lateness.quantile(0.9)), share_late_gt7d=float((lateness > 7).mean()),
            secondary_peaks_days_after=";".join(map(str, sec)) if sec else "",
            share_edited=float(g["edited"].mean()),
        ))
    return pd.DataFrame(rows)


def hazard_by_day(g: pd.DataFrame, s_max: int = 30, s_min: int = -10) -> pd.DataFrame:
    """Riesgo empírico por día en tiempo inverso: h(s) = n(entregados en el día s
    antes del plazo) / n(no entregados aún al empezar ese día), s = 30..-10."""
    s = np.floor(g["t_before"]).astype(int)          # día antes del plazo (0 = últimas 24 h; -1 = primer día tarde)
    total = len(g)
    out = []
    for k in range(s_max, s_min - 1, -1):
        at_risk = int((s <= k).sum())                # aún no entregados al comenzar el día k
        n_k = int((s == k).sum())
        out.append(dict(day_before=k, at_risk=at_risk, n=n_k, hazard=n_k / at_risk if at_risk else np.nan,
                        density=n_k / total))
    return pd.DataFrame(out)


def e1_figure(df: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    cmap = plt.get_cmap("viridis")
    for i, y in enumerate(YEARS):
        g = df[df["year"] == y]
        h = hazard_by_day(g, 25, -8)
        col = cmap(i / (len(YEARS) - 1))
        lab = f"{y} (T={C.window_days(y, True):.0f})" if C.window_days(y) else f"{y}"
        axes[0].plot(h["day_before"], 100 * h["density"], color=col, lw=1.4, label=lab)
        axes[1].plot(h["day_before"], 100 * h["hazard"], color=col, lw=1.4)
    for ax in axes:
        ax.invert_xaxis()
        ax.axvspan(S0, 0, color="grey", alpha=0.15, lw=0)
        ax.axvline(0, color="k", lw=0.8, ls="--")
        ax.set_xlabel("Days before the effective deadline (negative = late)")
    axes[0].set_ylabel("Share of reviews submitted that day (%)")
    axes[0].set_title("(a) Delivery distribution")
    axes[1].set_ylabel("Daily hazard (% of not-yet-submitted reviews)")
    axes[1].set_title("(b) Empirical hazard")
    axes[0].legend(fontsize=7, ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(C.FIG_DIR / "fig1_delivery_hazard.pdf")
    fig.savefig(C.FIG_DIR / "fig1_delivery_hazard.png", dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------- E2 ---------
def fit_trunc_exp(u: np.ndarray, U: float) -> float:
    """MLE de la tasa de una exponencial truncada en (0, U]: media = 1/g - U/(e^{gU}-1)."""
    m = float(np.mean(u))
    f = lambda g: 1.0 / g - U / np.expm1(g * U) - m
    with np.errstate(over="ignore"):
        return brentq(f, 1e-4, 20.0)


def loglik_trunc_exp(u: np.ndarray, g: float, U: float) -> float:
    return float(np.sum(np.log(g) - g * u - np.log(-np.expm1(-g * U))))


def fit_mixture(g: pd.DataFrame, U: float | None = None) -> dict:
    """Mezcla: pre-plazo anclado (exp. inversa en s-S0, truncada en U), átomo del día
    del plazo, cola tardía exponencial. Devuelve parámetros y verosimilitudes."""
    s = g["t_before"].values
    pre = s[s > S0]
    late = -s[s < 0]
    n = len(s)
    if U is None:
        U = float(np.quantile(pre, 0.995)) - S0
    u = np.clip(pre - S0, 1e-6, U)
    gamma = fit_trunc_exp(u, U)
    mu = 1.0 / max(late.mean(), 1e-6)
    out = dict(n=n, share_pre=len(pre) / n, share_atom=float(((s >= 0) & (s <= S0)).mean()), share_late=len(late) / n,
               gamma=gamma, halflife_days=np.log(2) / gamma, mean_days_before_pre=float(pre.mean()),
               mu=mu, mean_lateness=float(late.mean()), U=U)
    out["ll_anchored"] = loglik_trunc_exp(u, gamma, U)
    # modelo alternativo (P1 original): riesgo constante desde la asignación => densidad
    # decreciente en x = U - u (misma muestra pre-plazo, mismo soporte)
    x = U - u
    if x.mean() < U / 2:
        lam = fit_trunc_exp(x, U)
        out["ll_const"] = loglik_trunc_exp(x, lam, U)
    else:  # la densidad crece hacia el plazo: el MLE con riesgo constante (lambda>=0) es el límite uniforme
        lam = 0.0
        out["ll_const"] = float(-len(x) * np.log(U))
    out["lambda_const"] = lam
    out["llr_anchored_vs_const"] = out["ll_anchored"] - out["ll_const"]
    out["llr_per_obs"] = out["llr_anchored_vs_const"] / len(u)
    return out


def e2(df: pd.DataFrame, n_boot: int = 300) -> pd.DataFrame:
    rows = []
    for y in YEARS:
        g = df[df["year"] == y]
        U = (C.window_days(y, True) - S0) if C.window_days(y) else None
        base = fit_mixture(g, U)
        keys = ["share_pre", "share_atom", "share_late", "gamma", "halflife_days", "mu", "mean_lateness"]
        est, lo, hi = cluster_bootstrap(g, lambda d: [fit_mixture(d, U)[k] for k in keys], n_boot=n_boot, seed=C.SEED + y)
        row = dict(year=y, window=C.window_days(y, True), n=base["n"], U=base["U"])
        for k, e, l, h in zip(keys, est, lo, hi):
            row[k] = e; row[k + "_lo"] = l; row[k + "_hi"] = h
        row["ll_anchored"] = base["ll_anchored"]; row["ll_const"] = base["ll_const"]
        row["llr_per_obs"] = base["llr_per_obs"]; row["lambda_const"] = base["lambda_const"]
        rows.append(row)
        print(f"E2 {y}: pre {row['share_pre']:.3f} atom {row['share_atom']:.3f} late {row['share_late']:.3f} "
              f"gamma {row['gamma']:.3f} [{row['gamma_lo']:.3f},{row['gamma_hi']:.3f}] (half-life {row['halflife_days']:.2f} d) "
              f"mu {row['mu']:.3f}; LLR/obs anchored vs const = {row['llr_per_obs']:.3f}")
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- E3 ---------
QS = [0.25, 0.5, 0.75, 0.9]


def e3_stats(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for y in YEARS:
        T = C.window_days(y, True)
        if T is None:
            continue
        g = df[(df["year"] == y)]
        ontime = g[~g["late"]]
        r = dict(year=y, T=T, n=len(g))
        for q in QS:
            r[f"s_q{int(q*100)}"] = float(ontime["t_before"].quantile(q))          # días antes del plazo
            r[f"x_q{int(q*100)}"] = float(ontime["t_since_assign"].quantile(1 - q))  # días desde la asignación
        r["share_atom"] = float((g["phase"] == "atom").mean())
        r["share_late"] = float(g["late"].mean())
        r["share_first_week"] = float((ontime["t_since_assign"] <= 7).mean())
        r["gamma"] = fit_mixture(g, T - S0)["gamma"]
        rows.append(r)
    return pd.DataFrame(rows)


def slope(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    return float(np.polyfit(x, y, 1)[0])


def e3(df: pd.DataFrame, n_boot: int = 300) -> tuple[pd.DataFrame, pd.DataFrame]:
    tab = e3_stats(df)
    cols = [c for c in tab.columns if c.startswith(("s_q", "x_q"))] + ["share_atom", "share_late", "gamma"]
    # bootstrap por manuscrito dentro de cada año para la pendiente sobre T
    g = rng(C.SEED + 3)
    sub = df[df["t_since_assign"].notna()]
    groups = {y: {k: np.asarray(v) for k, v in d.groupby("forum").groups.items()} for y, d in sub.groupby("year")}
    reps = []
    for _ in range(n_boot):
        idx = []
        for y, gr in groups.items():
            ids = list(gr.keys())
            pick = g.choice(len(ids), size=len(ids), replace=True)
            idx.append(np.concatenate([gr[ids[i]] for i in pick]))
        t = e3_stats(sub.loc[np.concatenate(idx)])
        reps.append([slope(t["T"], t[c]) for c in cols])
    reps = np.array(reps)
    slopes = pd.DataFrame({"stat": cols, "slope_per_day": [slope(tab["T"], tab[c]) for c in cols],
                           "lo": np.percentile(reps, 2.5, axis=0), "hi": np.percentile(reps, 97.5, axis=0)})
    slopes["parkinson_null"] = [0.0 if c.startswith("s_q") else (1.0 if c.startswith("x_q") else np.nan) for c in cols]
    # error estándar clásico de la regresión entre años (n = 7 ventanas): la incertidumbre
    # entre cohortes no la captura el bootstrap por manuscrito
    ses = []
    for c in cols:
        X = np.column_stack([np.ones(len(tab)), tab["T"].values]); yv = tab[c].values
        beta, res, *_ = np.linalg.lstsq(X, yv, rcond=None)
        resid = yv - X @ beta
        s2 = (resid @ resid) / (len(tab) - 2)
        ses.append(float(np.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1])))
    slopes["se_between_years"] = ses
    return tab, slopes


def e3_figure(tab: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for q, m in zip(QS, ["o", "s", "^", "d"]):
        axes[0].plot(tab["T"], tab[f"s_q{int(q*100)}"], marker=m, ls="", label=f"q{int(q*100)}")
        axes[1].plot(tab["T"], tab[f"x_q{int(q*100)}"], marker=m, ls="", label=f"q{int(100-q*100)}")
    for _, r in tab.iterrows():
        axes[0].annotate(str(int(r["year"])), (r["T"], r["s_q50"]), fontsize=7, xytext=(3, 3), textcoords="offset points")
    xs = np.array([11, 23])
    axes[1].plot(xs, xs - 1, color="grey", ls="--", lw=0.8, label="T − 1 (full Parkinson)")
    axes[0].set_xlabel("Effective review window T (days)"); axes[0].set_ylabel("Days before the deadline (on-time reviews)")
    axes[0].set_title("(a) Quantiles of days before the deadline vs T")
    axes[1].set_xlabel("Effective review window T (days)"); axes[1].set_ylabel("Days since assignment (on-time reviews)")
    axes[1].set_title("(b) Quantiles of days since assignment vs T")
    axes[0].legend(fontsize=8, frameon=False); axes[1].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(C.FIG_DIR / "fig2_parkinson.pdf"); fig.savefig(C.FIG_DIR / "fig2_parkinson.png", dpi=200)
    plt.close(fig)


def main(n_boot: int = 300):
    df = load_reviews()
    t1 = e1_table(df)
    to_latex(t1, C.TAB_DIR / "table1_delivery", "Delivery of ICLR reviews relative to the deadline, 2018–2026", "tab:delivery", "%.3f")
    print(t1[["year", "window_official", "n_reviews", "median_days_before", "share_pre", "share_deadline_date", "share_late", "mean_lateness_days", "share_late_gt7d", "secondary_peaks_days_after"]].round(3).to_string())
    e1_figure(df)
    haz = pd.concat([hazard_by_day(g, 30, -10).assign(year=y) for y, g in df.groupby("year")])
    haz.to_csv(C.TAB_DIR / "hazard_by_day.csv", index=False)
    save_result("E1", t1.to_dict(orient="records"))

    t2 = e2(df, n_boot)
    to_latex(t2, C.TAB_DIR / "table2_mixture", "Mixture estimates by year", "tab:mixture", "%.3f")
    save_result("E2", t2.to_dict(orient="records"))

    t3, slopes = e3(df, n_boot)
    to_latex(t3, C.TAB_DIR / "table_e3_quantiles", "Quantiles of delivery time vs review window", "tab:parkinson", "%.2f")
    to_latex(slopes, C.TAB_DIR / "table_e3_slopes", "Slopes of delivery quantiles on the window T", "tab:parkinson_slopes", "%.3f")
    print(t3.round(2).to_string()); print(slopes.round(3).to_string())
    e3_figure(t3)
    save_result("E3", dict(table=t3.to_dict(orient="records"), slopes=slopes.to_dict(orient="records")))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 300)
