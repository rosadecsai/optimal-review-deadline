"""E4 calidad del informe por tramo temporal (efectos fijos de manuscrito);
E5 compresión entre ventanas anuales (12–22 días); E6 tomografía del consenso
por tramo temporal (kappa por pares de tramos, (f, s, c) con paneles triádicos y
tomografía diferencial de los informes de última hora y tardíos)."""
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
from utils import cluster_bootstrap, fe_regression, load_reviews, save_result, to_latex, rng  # noqa: E402

YEARS = sorted(C.YEARS)
BINS = C.BIN_ORDER            # early, mid, last48h, last24h, late
OUTCOMES = {
    "log_words": "log words (all text fields)",
    "bullets": "bullet/numbered items",
    "qmarks": "question marks",
    "confidence": "self-reported confidence (1–5)",
    "score": "overall recommendation",
    "dev_panel": "|score − mean of co-reviewers|",
    "words_weak": "words in weaknesses (2023–2026)",
    "words_quest": "words in questions (2024–2026)",
    "thoroughness": "thoroughness in reading (ICLR 2020, 1–3)",
    "check_deriv": "checked derivations (ICLR 2020, 0–2)",
    "check_exper": "checked experiments (ICLR 2020, 0–2)",
}


def add_bin_dummies(df: pd.DataFrame) -> list[str]:
    cols = []
    for b in BINS[1:]:
        df[f"bin_{b}"] = (df["bin"] == b).astype(float)
        cols.append(f"bin_{b}")
    return cols


# ---------------------------------------------------------------- E4 ---------
def e4(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    d = df.copy()
    cols = add_bin_dummies(d)
    rows, byyear = [], []
    for y_, lab in OUTCOMES.items():
        sub = d[d[y_].notna()]
        if sub.empty:
            continue
        res = fe_regression(sub, y_, cols)
        if res is None:
            continue
        for b in BINS[1:]:
            k = f"bin_{b}"
            if k in res.index:
                rows.append(dict(outcome=y_, label=lab, bin=b, coef=res.loc[k, "coef"], se=res.loc[k, "se"], p=res.loc[k, "p"],
                                 n=res.attrs["n"], n_papers=res.attrs["n_groups"], mean_early=float(sub.loc[sub["bin"] == "early", y_].mean())))
        for yr, g in sub.groupby("year"):
            r2 = fe_regression(g, y_, cols)
            if r2 is None:
                continue
            for b in BINS[1:]:
                k = f"bin_{b}"
                if k in r2.index:
                    byyear.append(dict(outcome=y_, year=yr, bin=b, coef=r2.loc[k, "coef"], se=r2.loc[k, "se"], n=r2.attrs["n"]))
    tab = pd.DataFrame(rows)
    tab["pct_effect"] = np.where(tab["outcome"] == "log_words", 100 * (np.exp(tab["coef"]) - 1), np.nan)
    return tab, pd.DataFrame(byyear)


def e4_figure(df: pd.DataFrame, byyear: pd.DataFrame):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, (out, ttl) in zip(axes, [("log_words", "(a) Length: log words vs early reviews"),
                                     ("confidence", "(b) Confidence vs early reviews"),
                                     ("dev_panel", "(c) Disagreement with co-reviewers vs early")]):
        t = byyear[byyear["outcome"] == out]
        for b, m in zip(BINS[1:], ["o", "s", "^", "x"]):
            tb = t[t["bin"] == b].sort_values("year")
            ax.errorbar(tb["year"], tb["coef"], yerr=1.96 * tb["se"], marker=m, ms=4, lw=1, capsize=2, label=b)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_title(ttl, fontsize=10); ax.set_xlabel("ICLR year")
        ax.set_ylabel("Within-manuscript coefficient")
    axes[0].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(C.FIG_DIR / "fig3_quality_by_bin.pdf"); fig.savefig(C.FIG_DIR / "fig3_quality_by_bin.png", dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------- E5 ---------
def year_stats(g: pd.DataFrame) -> list[float]:
    early = g[g["bin"] == "early"]
    last = g[g["bin"] == "last24h"]
    return [g["words"].mean(), g["confidence"].mean(), g["dev_panel"].mean(), g["score"].std(),
            (g["phase"] == "atom").mean(), g["late"].mean(), (-g.loc[g["late"], "t_before"]).mean(),
            np.log(last["words"].clip(lower=1)).mean() - np.log(early["words"].clip(lower=1)).mean(),
            last["confidence"].mean() - early["confidence"].mean(),
            last["dev_panel"].mean() - early["dev_panel"].mean()]


E5_KEYS = ["words", "confidence", "dev_panel", "sd_score", "share_atom", "share_late", "mean_lateness",
           "penalty_logwords_last24h", "penalty_conf_last24h", "penalty_dev_last24h"]


def e5(df: pd.DataFrame, n_boot: int = 300) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for y in YEARS:
        g = df[df["year"] == y]
        est, lo, hi = cluster_bootstrap(g, year_stats, n_boot=n_boot, seed=C.SEED + 5 + y)
        r = dict(year=y, T=C.window_days(y), T_eff=C.window_days(y, True), n=len(g), form="2018-19" if y <= 2019 else ("2020" if y == 2020 else ("2021" if y == 2021 else ("2022-23" if y <= 2023 else "2024-26"))))
        for k, e, l, h in zip(E5_KEYS, est, lo, hi):
            r[k] = e; r[k + "_lo"] = l; r[k + "_hi"] = h
        rows.append(r)
    tab = pd.DataFrame(rows)
    # contrastes: 2023 (12 d) frente a 2022 (misma plantilla, 18/19 d); 2024-26 (20-21 d) como placebo
    def diff(a, b, k):
        ra, rb = tab[tab["year"] == a].iloc[0], tab[tab["year"] == b].iloc[0]
        se = np.sqrt(((ra[k + "_hi"] - ra[k + "_lo"]) / 3.92) ** 2 + ((rb[k + "_hi"] - rb[k + "_lo"]) / 3.92) ** 2)
        return ra[k] - rb[k], se
    contrasts = []
    for k in E5_KEYS:
        for a, b in [(2023, 2022), (2023, 2024), (2025, 2024), (2026, 2025), (2021, 2020)]:
            dlt, se = diff(a, b, k)
            contrasts.append(dict(stat=k, years=f"{a} vs {b}", dT=C.window_days(a, True) - C.window_days(b, True), diff=dlt, se=se, z=dlt / se if se > 0 else np.nan))
    return tab, pd.DataFrame(contrasts)


# ---------------------------------------------------------------- E6 ---------
ACCEPT_THRESHOLD = 6.0


def panel_moments(votes: pd.DataFrame) -> tuple[float, float, float, int]:
    """Momentos m1, m2, m3 con peso igual por panel; votes: columnas forum, accept."""
    g = votes.groupby("forum")["accept"].agg(["sum", "size"])
    K, n = g["sum"].values.astype(float), g["size"].values.astype(float)
    m1 = float(np.mean(K / n))
    ok2 = n >= 2
    m2 = float(np.mean((K[ok2] * (K[ok2] - 1)) / (n[ok2] * (n[ok2] - 1))))
    ok3 = n >= 3
    m3 = float(np.mean((K[ok3] * (K[ok3] - 1) * (K[ok3] - 2)) / (n[ok3] * (n[ok3] - 1) * (n[ok3] - 2)))) if ok3.any() else np.nan
    return m1, m2, m3, int(ok3.sum())


def invert_three_moments(m1: float, m2: float, m3: float) -> dict:
    """Tomografía del consenso (dos átomos): S=(m3−m1m2)/(m2−m1²), P=(m1m3−m2²)/(m2−m1²);
    átomos raíces de z²−Sz+P (a>b); c=a, s=1−b/a, f=(a−m1)/(a−b); kappa=(m2−m1²)/(m1(1−m1))."""
    var = m2 - m1 ** 2
    kappa = var / (m1 * (1 - m1))
    if not np.isfinite(m3) or var <= 0:
        return dict(kappa=kappa, f=np.nan, s=np.nan, c=np.nan, a=np.nan, b=np.nan)
    S = (m3 - m1 * m2) / var
    P = (m1 * m3 - m2 ** 2) / var
    disc = S ** 2 - 4 * P
    if disc < 0:
        return dict(kappa=kappa, f=np.nan, s=np.nan, c=np.nan, a=np.nan, b=np.nan)
    a = (S + np.sqrt(disc)) / 2; b = (S - np.sqrt(disc)) / 2
    f = (a - m1) / (a - b) if a != b else np.nan
    return dict(kappa=kappa, f=f, s=1 - b / a if a else np.nan, c=a, a=a, b=b)


def pair_kappa(v1: pd.DataFrame, v2: pd.DataFrame, same: bool) -> tuple[float, int]:
    """Kappa entre dos grupos de informes del mismo manuscrito (pares cruzados, o
    pares dentro del grupo si same=True)."""
    if same:
        g = v1.groupby("forum")["accept"].agg(["sum", "size"]); g = g[g["size"] >= 2]
        K, n = g["sum"].values.astype(float), g["size"].values.astype(float)
        pairs = n * (n - 1) / 2
        both = K * (K - 1) / 2; neither = (n - K) * (n - K - 1) / 2
        po = float((both + neither).sum() / pairs.sum()); p1 = p2 = float((K * (n - 1)).sum() / (n * (n - 1)).sum())
        npairs = int(pairs.sum())
    else:
        a = v1.groupby("forum")["accept"].agg(["sum", "size"]); b = v2.groupby("forum")["accept"].agg(["sum", "size"])
        j = a.join(b, lsuffix="1", rsuffix="2", how="inner")
        K1, n1, K2, n2 = [j[c].values.astype(float) for c in ["sum1", "size1", "sum2", "size2"]]
        pairs = n1 * n2
        agree = K1 * K2 + (n1 - K1) * (n2 - K2)
        po = float(agree.sum() / pairs.sum())
        p1 = float((K1 * n2).sum() / pairs.sum()); p2 = float((K2 * n1).sum() / pairs.sum())
        npairs = int(pairs.sum())
    pe = p1 * p2 + (1 - p1) * (1 - p2)
    return (po - pe) / (1 - pe), npairs


def differential_tomography(votes: pd.DataFrame, ref: str = "atom") -> dict:
    """(f, s, c) del grupo de referencia con sus tríadas; para los demás grupos
    (a_g, b_g) desde m1_g y el momento cruzado m2(g, ref) con f común:
      m1_g = f b_g + (1−f) a_g ;  m2(g,ref) = f b_g b_ref + (1−f) a_g a_ref."""
    out = {}
    vr = votes[votes["grp"] == ref]
    m1, m2, m3, n3 = panel_moments(vr)
    base = invert_three_moments(m1, m2, m3)
    out[ref] = dict(base, m1=m1, m2=m2, m3=m3, n_panels3=n3, method="triads")
    f, aE, bE = base["f"], base["a"], base["b"]
    A = vr.groupby("forum")["accept"].agg(["sum", "size"])
    for grp in [g_ for g_ in ["early", "mid", "atom", "late"] if g_ != ref]:
        vg = votes[votes["grp"] == grp]
        if vg.empty:
            continue
        m1g = float(vg.groupby("forum")["accept"].mean().mean())
        B = vg.groupby("forum")["accept"].agg(["sum", "size"])
        j = A.join(B, lsuffix="E", rsuffix="G", how="inner")
        m2gE = float(np.mean((j["sumE"] / j["sizeE"]) * (j["sumG"] / j["sizeG"])))
        rec = dict(m1=m1g, m2_cross=m2gE, n_panels_cross=len(j), method="differential")
        if np.isfinite(f) and np.isfinite(aE) and aE != bE:
            M = np.array([[1 - f, f], [(1 - f) * aE, f * bE]])
            ag, bg = np.linalg.solve(M, np.array([m1g, m2gE]))
            rec.update(a=ag, b=bg, c=ag, s=1 - bg / ag if ag else np.nan, f=f)
        # tríadas propias del grupo (f_g libre): comprobación de selección
        m1s, m2s, m3s, n3s = panel_moments(vg)
        own = invert_three_moments(m1s, m2s, m3s)
        rec.update({f"own_{k}": v for k, v in own.items()}, own_n_panels3=n3s)
        out[grp] = rec
    return out


MIN_PANELS3 = 500  # tríadas mínimas para invertir momentos con un grupo


def e6(df: pd.DataFrame, n_boot: int = 200) -> tuple[pd.DataFrame, pd.DataFrame]:
    d = df[df["score"].notna()].copy()
    d["accept"] = (d["score"] >= ACCEPT_THRESHOLD).astype(int)
    d["grp"] = d["bin"].map({"early": "early", "mid": "mid", "last48h": "atom", "last24h": "atom", "late": "late"})
    kap_rows, tomo_rows = [], []
    for y in YEARS:
        g = d[d["year"] == y]
        groups = {k: g[g["grp"] == k] for k in ["early", "mid", "atom", "late"]}
        kall, nall = pair_kappa(g, None, True)
        kap_rows.append(dict(year=y, pair="all", kappa=kall, n_pairs=nall, n_reviews=len(g), accept_rate=float(g["accept"].mean())))
        for i, a in enumerate(["early", "mid", "atom", "late"]):
            for b in ["early", "mid", "atom", "late"][i:]:
                k, n = pair_kappa(groups[a], groups[b], a == b)
                kap_rows.append(dict(year=y, pair=f"{a}-{b}", kappa=k, n_pairs=n))
        # kappa por pares con IC bootstrap (todas, early-early, atom-atom, early-atom, late-late, early-late)
        pairs = [("early", "early"), ("atom", "atom"), ("late", "late"), ("early", "atom"), ("early", "late"), ("atom", "late")]
        def kstat(dd):
            gg = {k: dd[dd["grp"] == k] for k in ["early", "mid", "atom", "late"]}
            return [pair_kappa(dd, None, True)[0]] + [pair_kappa(gg[a], gg[b], a == b)[0] for a, b in pairs]
        est, lo, hi = cluster_bootstrap(g, kstat, n_boot=n_boot, seed=C.SEED + 6 + y)
        row = dict(year=y, n_reviews=len(g), accept_rate=float(g["accept"].mean()))
        for name, e, l, h in zip(["all"] + [f"{a}-{b}" for a, b in pairs], est, lo, hi):
            row[f"kappa_{name}"] = e; row[f"kappa_{name}_lo"] = l; row[f"kappa_{name}_hi"] = h
        # tomografía por grupo con tríadas propias (solo grupos con >= MIN_PANELS3 tríadas)
        for grp in ["early", "atom", "late", "all"]:
            vg = g if grp == "all" else groups[grp]
            m1, m2, m3, n3 = panel_moments(vg)
            if n3 < MIN_PANELS3:
                continue
            def tstat(dd):
                mm = panel_moments(dd)
                t = invert_three_moments(*mm[:3])
                return [t["f"], t["s"], t["c"], t["kappa"]]
            te, tl, th = cluster_bootstrap(vg, tstat, n_boot=n_boot, seed=C.SEED + 60 + y)
            row[f"n3_{grp}"] = n3
            for k, e, l, h in zip(["f", "s", "c", "kappa"], te, tl, th):
                row[f"{grp}_{k}"] = e; row[f"{grp}_{k}_lo"] = l; row[f"{grp}_{k}_hi"] = h
        tomo_rows.append(row)
        msg = " ".join(f"{grp}:(f={row.get(grp+'_f', np.nan):.2f},s={row.get(grp+'_s', np.nan):.2f},c={row.get(grp+'_c', np.nan):.2f},n3={row.get('n3_'+grp, 0)})" for grp in ["all", "early", "atom", "late"] if f"n3_{grp}" in row)
        print(f"E6 {y}: kappa all {row['kappa_all']:.3f} early-early {row['kappa_early-early']:.3f} atom-atom {row['kappa_atom-atom']:.3f} early-late {row['kappa_early-late']:.3f} | {msg}")
    return pd.DataFrame(kap_rows), pd.DataFrame(tomo_rows)


def main(n_boot: int = 300):
    df = load_reviews()
    df = df[df["n_reviews"] >= 2]
    t4, by4 = e4(df)
    to_latex(t4, C.TAB_DIR / "table3_quality_by_bin", "Review quality proxies by submission timing (manuscript fixed effects)", "tab:quality", "%.3f")
    by4.to_csv(C.TAB_DIR / "table3_quality_by_bin_byyear.csv", index=False)
    print(t4[["outcome", "bin", "coef", "se", "p", "pct_effect", "n", "mean_early"]].round(3).to_string())
    e4_figure(df, by4)
    save_result("E4", dict(pooled=t4.to_dict(orient="records"), by_year=by4.to_dict(orient="records")))

    t5, c5 = e5(df, n_boot)
    to_latex(t5, C.TAB_DIR / "table4_windows", "Review outcomes by cohort and review window", "tab:windows", "%.3f")
    to_latex(c5, C.TAB_DIR / "table4_contrasts", "Contrasts between cohorts", "tab:contrasts", "%.3f")
    print(t5[["year", "T", "words", "confidence", "dev_panel", "share_atom", "share_late", "mean_lateness", "penalty_logwords_last24h", "penalty_conf_last24h", "penalty_dev_last24h"]].round(3).to_string())
    print(c5.round(3).to_string())
    save_result("E5", dict(table=t5.to_dict(orient="records"), contrasts=c5.to_dict(orient="records")))

    k6, t6 = e6(df, max(100, n_boot // 2))
    to_latex(k6, C.TAB_DIR / "table5_kappa_pairs", "Agreement (kappa) between reviews by timing", "tab:kappa_pairs", "%.3f")
    to_latex(t6, C.TAB_DIR / "table5_tomography", "Consensus tomography by timing", "tab:tomography", "%.3f")
    print(k6.pivot(index="year", columns="pair", values="kappa").round(3).to_string())
    save_result("E6", dict(kappa=k6.to_dict(orient="records"), tomography=t6.to_dict(orient="records")))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 300)
