"""Genera las tablas LaTeX del manuscrito (manuscript/tables/*.tex) a partir de
output/results.json y output/tables/*.csv, para que las cifras del texto y de las
tablas salgan del mismo sitio."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402

OUT = C.ROOT / "manuscript" / "tables"
OUT.mkdir(parents=True, exist_ok=True)
R = json.loads((C.ROOT / "output" / "results.json").read_text())
T = C.TAB_DIR

FORMS = {
    2018: ("review, rating 1--10, confidence", "v1"), 2019: ("review, rating 1--10, confidence", "v1"),
    2020: ("review, rating 1/3/6/8, self-reported effort", "v1"), 2021: ("review, rating 1--10, confidence", "v1"),
    2022: ("main review, recommendation, correctness, novelty", "v1"),
    2023: ("strengths/weaknesses, recommendation, correctness, novelty", "v1"),
    2024: ("summary, strengths, weaknesses, questions, rating, soundness", "v2"),
    2025: ("idem, numeric confidence", "v2"), 2026: ("idem, rating 0--10 in steps of 2", "v2"),
}


def fmt(x, d=2):
    return "--" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.{d}f}"


def pct(x, d=0):
    return f"{100 * x:.{d}f}"


def wrap(body: str, caption: str, label: str, cols: str, note: str = "", size: str = "\\small", env: str = "table", colsep: str = "4pt", width: str = "\\textwidth") -> str:
    s = [f"\\begin{{{env}}}[{'p' if env == 'sidewaystable' else 't'}]", f"\\caption{{{caption}}}\\label{{{label}}}", size, f"\\setlength{{\\tabcolsep}}{{{colsep}}}",
         f"\\begin{{tabular*}}{{{width}}}{{@{{\\extracolsep{{\\fill}}}}{cols}@{{}}}}", "\\toprule", body, "\\bottomrule", "\\end{tabular*}"]
    if note:
        s.append(f"\\footnotetext{{{note}}}")
    s.append(f"\\end{{{env}}}")
    return "\n".join(s) + "\n"


# ---------------------------------------------------------------- Table 1 ---
def table1():
    e1 = pd.DataFrame(R["E1"]).set_index("year")
    rows = []
    for y in sorted(C.YEARS):
        a = C.YEARS[y]["assign"]; d = C.YEARS[y]["deadline"]; e = C.EFFECTIVE_DEADLINE[y]
        dl = e.strftime("%d %b %Y") + ("" if e == d else "$^{\\dagger}$")
        rows.append(" & ".join([
            str(y), a.strftime("%d %b %Y") if a else "n/a", dl,
            fmt(C.window_days(y, True), 0) if C.window_days(y) else "--",
            f"{int(e1.loc[y, 'n_papers']):,}", f"{int(e1.loc[y, 'n_reviews']):,}", FORMS[y][0]]) + " \\\\")
    body = "Cohort & Assignment & Deadline (AoE) & $T$ & Submissions & Reviews & Review form \\\\\n\\midrule\n" + "\n".join(rows)
    note = ("Assignment dates and deadlines from the ICLR reviewer guides (iclr.cc); no assignment date is published for 2018--2019. "
            "$^{\\dagger}$Effective deadline, the day (Anywhere on Earth) after which the hourly submission rate collapses, where it differs from the published date: "
            "2022 (published 1 Nov 2021; one-day extension not recorded on the official page) and 2026 (published 1 Nov 2025; the invitation due date, 1 Nov 2025 11:59 UTC, is 31 Oct AoE). "
            "$T$ = days from assignment to the effective deadline. Submissions = manuscripts with at least one official review.")
    (OUT / "tab1_cohorts.tex").write_text(wrap(body, "The nine ICLR cohorts: review windows, sample sizes and review forms", "tab:cohorts", "llllrrp{3.1cm}", note, "\\scriptsize", colsep="3pt"))


# ---------------------------------------------------------------- Table 2 ---
def table2():
    e1 = pd.DataFrame(R["E1"]).set_index("year"); e2 = pd.DataFrame(R["E2"]).set_index("year")
    rows = []
    for y in sorted(C.YEARS):
        a = e1.loc[y]
        rows.append(" & ".join([str(y), fmt(a["median_days_before"], 2), pct(a["share_pre"]), pct(a["share_deadline_date"]), pct(a["share_last24h"]),
                                pct(a["share_late"]), fmt(a["mean_lateness_days"], 1), fmt(a["median_lateness_days"], 1), pct(a["share_late_gt7d"])]) + " \\\\")
    body = ("Cohort & Median days before & $>$48 h (\\%) & Last 48 h (\\%) & Last 24 h (\\%) & Late (\\%) & Mean lateness (d) & Median lateness (d) & Late $>$7 d (\\%) \\\\\n\\midrule\n" + "\n".join(rows))
    note = "Timing relative to the effective deadline (23:59 AoE). Shares over all reviews of the cohort; lateness statistics over late reviews."
    (OUT / "tab2_delivery.tex").write_text(wrap(body, "When reviews arrive: delivery relative to the deadline by cohort", "tab:delivery", "lrrrrrrrr", note, "\\footnotesize", colsep="3pt", env="sidewaystable", width="0.95\\textheight"))
    rows = []
    for y in sorted(C.YEARS):
        b = e2.loc[y]
        rows.append(" & ".join([str(y), pct(b["share_atom"], 1), pct(b["share_late"], 1),
                                f"{fmt(b['gamma'], 3)} ({fmt(b['gamma_lo'], 3)}--{fmt(b['gamma_hi'], 3)})", fmt(b["halflife_days"], 1),
                                f"{fmt(b['mu'], 3)} ({fmt(b['mu_lo'], 3)}--{fmt(b['mu_hi'], 3)})", fmt(b["mean_lateness"], 1), fmt(b["llr_per_obs"], 2)]) + " \\\\")
    body = ("Cohort & Atom, last 48 h (\\%) & Late tail (\\%) & $\\hat\\gamma$ per day (95\\% CI) & Half-life (d) & $\\hat\\mu$ per day (95\\% CI) & Mean lateness (d) & LLR per review \\\\\n\\midrule\n" + "\n".join(rows))
    note = ("Maximum-likelihood estimates of the three-component mixture: deadline-anchored component (reverse exponential with rate $\\gamma$ in days before the final 48 h, truncated at the window; half-life $=\\ln 2/\\hat\\gamma$), atom on the final 48 h, geometric late tail with daily rate $\\mu$. "
            "LLR per review: log-likelihood ratio per pre-deadline review of the deadline-anchored component against the constant-hazard-from-assignment alternative, whose maximum-likelihood limit is the uniform density in all nine cohorts. "
            "95\\% intervals from 500 cluster (manuscript) bootstrap replications.")
    (OUT / "tab2b_mixture.tex").write_text(wrap(body, "Mixture estimates: anchoring to the deadline, atom and late tail by cohort", "tab:mixture", "lrrlrlrr", note, "\\footnotesize", colsep="3pt", env="sidewaystable", width="0.95\\textheight"))


# ---------------------------------------------------------------- Table 3 ---
def table3():
    s = pd.DataFrame(R["E3"]["slopes"])
    labels = {"s_q25": "$Q_{25}$ of days before the deadline", "s_q50": "Median days before the deadline", "s_q75": "$Q_{75}$ of days before the deadline",
              "s_q90": "$Q_{90}$ of days before the deadline", "x_q25": "$Q_{75}$ of days since assignment", "x_q50": "Median days since assignment",
              "x_q75": "$Q_{25}$ of days since assignment", "x_q90": "$Q_{10}$ of days since assignment", "share_atom": "Share in the last 48 h",
              "share_late": "Share late", "gamma": "Anchoring rate $\\hat\\gamma$"}
    rows = []
    for _, r in s.iterrows():
        null = "--" if np.isnan(r["parkinson_null"]) else fmt(r["parkinson_null"], 0)
        rows.append(" & ".join([labels[r["stat"]], fmt(r["slope_per_day"], 3), f"[{fmt(r['lo'], 3)}, {fmt(r['hi'], 3)}]", fmt(r["se_between_years"], 3), null]) + " \\\\")
    body = "Statistic (on-time reviews) & Slope on $T$ & 95\\% CI & SE between cohorts & Parkinson value \\\\\n\\midrule\n" + "\n".join(rows)
    note = ("OLS slope of the cohort statistic on the effective window $T$ across the seven cohorts with a published assignment date (2020--2026, $T=12$--21 days). "
            "Bootstrap: manuscripts resampled within each cohort (500 replications); the between-cohort standard error treats the seven cohorts as the units. "
            "Under Parkinson's law the days-before-deadline quantiles do not depend on $T$ (slope 0) and the days-since-assignment quantiles move one-for-one (slope 1).")
    (OUT / "tab3_parkinson.tex").write_text(wrap(body, "Parkinson's law across review windows: slopes of delivery quantiles on the window $T$", "tab:parkinson", "lrlrr", note, "\\footnotesize", colsep="3pt"))


# ---------------------------------------------------------------- Table 4 ---
def table4():
    t = pd.DataFrame(R["E4"]["pooled"])
    order = ["log_words", "bullets", "qmarks", "words_weak", "words_quest", "confidence", "score", "dev_panel"]
    labels = {"log_words": "log words (all text fields)", "bullets": "Bullet / numbered items", "qmarks": "Question marks", "words_weak": "Words in weaknesses (2023--2026)",
              "words_quest": "Words in questions (2024--2026)", "confidence": "Confidence (1--5)", "score": "Overall score", "dev_panel": "$|$score $-$ mean of co-reviewers$|$"}
    rows = []
    for o in order:
        tt = t[t["outcome"] == o].set_index("bin")
        cells = [labels[o], fmt(tt["mean_early"].iloc[0], 2)]
        for b in ["mid", "last48h", "last24h", "late"]:
            c, se, p = tt.loc[b, "coef"], tt.loc[b, "se"], tt.loc[b, "p"]
            star = "$^{***}$" if p < 0.001 else ("$^{**}$" if p < 0.01 else ("$^{*}$" if p < 0.05 else ""))
            cells.append(f"{fmt(c, 3)}{star} ({fmt(se, 3)})")
        rows.append(" & ".join(cells) + " \\\\")
    n = int(t["n"].max()); npap = int(t["n_papers"].max())
    body = "Outcome & Mean, early & 2--4 days before & Last 48 h & Last 24 h & Late \\\\\n\\midrule\n" + "\n".join(rows)
    note = (f"Within-manuscript regressions on timing-bin dummies (reference: more than 4 days before the deadline), manuscript fixed effects, standard errors clustered by manuscript in parentheses; "
            f"$N={n:,}$ reviews in {npap:,} manuscripts with at least two reviews, 2018--2026 pooled. $^{{***}}p<0.001$, $^{{**}}p<0.01$, $^{{*}}p<0.05$. "
            "For log words the last-24-hour coefficient corresponds to $100(e^{-0.085}-1)=-8.2\\%$. Confidence is unavailable for 2020.")
    (OUT / "tab4_quality.tex").write_text(wrap(body, "Review quality proxies by submission timing (manuscript fixed effects)", "tab:quality", "lrllll", note, "\\scriptsize", colsep="2.5pt", env="sidewaystable", width="0.95\\textheight"))


# ---------------------------------------------------------------- Table 5 ---
def table5():
    k = pd.DataFrame(R["E6"]["kappa"])
    piv = k.pivot(index="year", columns="pair", values="kappa")
    npairs = k.pivot(index="year", columns="pair", values="n_pairs")
    cols = ["all", "early-early", "early-atom", "atom-atom", "early-late", "atom-late", "late-late"]
    rows = []
    for y in range(2018, 2026):
        rows.append(" & ".join([str(y), f"{int(npairs.loc[y, 'all']):,}"] + [fmt(piv.loc[y, c], 2) for c in cols] + [f"{int(npairs.loc[y, 'early-early']):,}"]) + " \\\\")
    body = "Cohort & Pairs & All & Early--early & Early--last 48 h & Last 48 h--last 48 h & Early--late & Last 48 h--late & Late--late & Early--early pairs \\\\\n\\midrule\n" + "\n".join(rows)
    note = ("Cohen's $\\kappa$ on the binary accept vote (score $\\geq 6$) between reviews of the same manuscript, by the timing of the two reviews (early: more than 4 days before the deadline). "
            "2026 is omitted: its new 0--10 scale in steps of 2 (``4: marginally below \\ldots would not mind if accepted'') gives $\\kappa=0.11$ in every timing pair and is not comparable. "
            "Triadic consensus-tomography estimates by timing group are reported in Online Resource~2.")
    (OUT / "tab5_kappa.tex").write_text(wrap(body, "Agreement with co-reviewers does not depend on when the review was written", "tab:kappa", "lrcccccccr", note, "\\scriptsize", colsep="3pt", env="sidewaystable", width="0.95\\textheight"))


# ---------------------------------------------------------------- Table 6 ---
def table6():
    t = pd.DataFrame(R["E7"]["tstar"])
    rows = []
    for _, r in t.iterrows():
        name = ("Conference" if r["conference"] else f"Journal, $L={r['L']}$")
        rows.append(" & ".join([name, f"{r['delta']:.4f}".rstrip("0"), f"{r['eta']:.1f}", f"{int(r['T_star'])}", fmt(r["a_Tstar"], 2), pct(r["share_atom_Tstar"] + r["share_late_Tstar"])] + [fmt(r[f"loss_T{d}"], 1) for d in (7, 14, 21, 28, 42)]) + " \\\\")
    body = "Scenario & $\\delta$ & $\\eta$ & $T^*$ & $a(T^*)$ & Deadline mass (\\%) & \\multicolumn{5}{c}{Loss at $T$ = 7 / 14 / 21 / 28 / 42 days (\\% of $\\bar q$)} \\\\\n\\midrule\n" + "\n".join(rows)
    note = ("$\\delta$: delay cost per day as a fraction of the value of a complete report; $\\eta$: compression elasticity, $\\Delta(T)=\\Delta\\,(20.5/T)^{\\eta}$ with $\\Delta=0.085$ (last-24-hour length penalty, Table~\\ref{tab:quality}); "
            "$L$: loss when the invitation is declined; deadline mass: share of reviews in the final 48 h or late at $T^*$. Conference: reviewers pre-committed ($a\\equiv1$, $L=0$). Journal: $a(T)=\\bar a(1-e^{-\\lambda_a T})$ through the two arms of \\citet{chetty2014} "
            "($\\bar a=0.684$, $\\lambda_a=0.099$). Delivery moments for every $T$ from the 2024--2026 hazard profile (Sect.~\\ref{sec:calibration}). 95\\% bootstrap intervals for $T^*$ have width 0--1 day and are omitted.")
    (OUT / "tab6_tstar.tex").write_text(wrap(body, "The optimal review window by scenario", "tab:tstar", "lllrrrrrrrr", note, "\\scriptsize", colsep="3pt", env="sidewaystable", width="0.95\\textheight"))


# ---------------------------------------------------------------- OR2 -------
def or2_tables():
    # E5 contrasts
    c = pd.DataFrame(R["E5"]["contrasts"])
    keep = ["words", "confidence", "dev_panel", "share_atom", "share_late", "mean_lateness", "penalty_logwords_last24h", "penalty_conf_last24h"]
    lab = {"words": "Mean words", "confidence": "Mean confidence", "dev_panel": "Mean deviation from panel", "share_atom": "Share last 48 h", "share_late": "Share late",
           "mean_lateness": "Mean lateness (days)", "penalty_logwords_last24h": "Last-24 h penalty, log words", "penalty_conf_last24h": "Last-24 h penalty, confidence"}
    rows = []
    for k in keep:
        cc = c[c["stat"] == k]
        rows.append(" & ".join([lab[k]] + [f"{fmt(r['diff'], 3)} ({fmt(r['z'], 1)})" for _, r in cc.iterrows()]) + " \\\\")
    cc0 = c[c["stat"] == "words"]
    heads = [r["years"].replace(" vs ", "--") for _, r in cc0.iterrows()]
    heads2 = [f"$\\Delta T={int(r['dT'])}$" for _, r in cc0.iterrows()]
    body = "Statistic & " + " & ".join(heads) + " \\\\\n & " + " & ".join(heads2) + " \\\\\n\\midrule\n" + "\n".join(rows)
    note = "Differences between cohort means (z-statistics from cluster-bootstrap standard errors in parentheses). 2023 vs 2022 compares the 12-day window with the 19-day window under the same review form; 2025 vs 2024 and 2026 vs 2025 are placebo contrasts (same form, one-day change)."
    (OUT / "tabS1_contrasts.tex").write_text(wrap(body, "Cohort contrasts: compression across review windows (E5)", "tab:contrasts", "lccccc", note, "\\scriptsize", colsep="3pt", env="sidewaystable", width="0.95\\textheight"))
    # tomography
    t = pd.DataFrame(R["E6"]["tomography"])
    rows = []
    for _, r in t.iterrows():
        cells = [str(int(r["year"]))]
        for g in ["all", "early", "atom", "late"]:
            if f"n3_{g}" in r and not np.isnan(r.get(f"n3_{g}", np.nan)):
                cells.append(f"{fmt(r[g + '_f'], 2)} / {fmt(r[g + '_s'], 2)} / {fmt(r[g + '_c'], 2)} ({int(r['n3_' + g]):,})")
            else:
                cells.append("--")
        rows.append(" & ".join(cells) + " \\\\")
    body = "Cohort & All reviews & Early ($>$4 d) & Last 48 h & Late \\\\\n\\midrule\n" + "\n".join(rows)
    note = ("Plug-in estimates $\\hat f/\\hat s/\\hat c$ (share of flawed manuscripts, selectivity, curation threshold) from the first three panel moments of the accept vote, computed within each timing group "
            "(number of panels with at least three reviews of that group in parentheses; groups with fewer than 500 such panels omitted). The moment inversion is biased upwards in finite samples "
            "(percentile bootstrap intervals lie above the plug-in values), so only the qualitative comparison across groups is informative.")
    (OUT / "tabS2_tomography.tex").write_text(wrap(body, "Consensus tomography by timing group (E6)", "tab:tomography", "lllll", note, "\\scriptsize", colsep="3pt", env="sidewaystable", width="0.95\\textheight"))
    # E4 by year (log words and confidence, last 24 h and late)
    b = pd.DataFrame(R["E4"]["by_year"])
    rows = []
    for y in sorted(b["year"].unique()):
        cells = [str(int(y))]
        for o in ["log_words", "confidence", "score", "dev_panel"]:
            for bn in ["last24h", "late"]:
                s = b[(b["year"] == y) & (b["outcome"] == o) & (b["bin"] == bn)]
                cells.append("--" if s.empty else f"{fmt(s['coef'].iloc[0], 2)} ({fmt(s['se'].iloc[0], 2)})")
        rows.append(" & ".join(cells) + " \\\\")
    body = "Cohort & \\multicolumn{2}{c}{log words} & \\multicolumn{2}{c}{Confidence} & \\multicolumn{2}{c}{Score} & \\multicolumn{2}{c}{Deviation from panel} \\\\\n & Last 24 h & Late & Last 24 h & Late & Last 24 h & Late & Last 24 h & Late \\\\\n\\midrule\n" + "\n".join(rows)
    note = "Within-manuscript coefficients (reference: more than 4 days before the deadline), by cohort; clustered standard errors in parentheses."
    (OUT / "tabS3_quality_by_year.tex").write_text(wrap(body, "Timing penalties by cohort (E4)", "tab:quality_year", "lcccccccc", note, "\\scriptsize", colsep="2.5pt", env="sidewaystable", width="0.95\\textheight"))
    # validation
    v = pd.DataFrame(R["E7"]["validation"])
    rows = []
    for _, r in v.iterrows():
        rows.append(" & ".join([str(int(r["year"])), fmt(r["T"], 0), "yes" if r["in_calibration"] else "no",
                                f"{pct(r['obs_share_atom'])} / {pct(r['pred_share_atom'])}", f"{pct(r['obs_share_late'])} / {pct(r['pred_share_late'])}",
                                f"{pct(r['obs_share_pre'])} / {pct(r['pred_share_pre'])}", f"{fmt(r['obs_mean_lateness'], 1)} / {fmt(r['pred_mean_lateness'], 1)}"]) + " \\\\")
    body = "Cohort & $T$ & Calibration & Last 48 h (\\%) obs/pred & Late (\\%) obs/pred & $>$48 h (\\%) obs/pred & Lateness (d) obs/pred \\\\\n\\midrule\n" + "\n".join(rows)
    note = "Predictions for each cohort's window from the hazard profile estimated on the 2024--2026 cohorts (Sect.~6.1 of the paper)."
    (OUT / "tabS4_validation.tex").write_text(wrap(body, "Out-of-sample prediction of cohort timing from the 2024--2026 hazard profile", "tab:validation", "lrlcccc", note, "\\scriptsize", colsep="3pt", env="sidewaystable", width="0.95\\textheight"))
    # sensitivity
    s = pd.DataFrame(R["E7"]["sensitivity"])
    pj = s.pivot_table(index=["Delta", "eta"], columns=["lam_a", "delta"], values="T_star_journal")
    pc = s.pivot_table(index=["Delta", "eta"], columns="delta", values="T_star_conference")
    rows = []
    for (D, eta), r in pj.iterrows():
        rows.append(" & ".join([fmt(D, 3), fmt(eta, 1)] + [str(int(r[c])) for c in pj.columns] + [str(int(pc.loc[(D, eta), c])) for c in pc.columns]) + " \\\\")
    lam_labels = [f"$\\lambda_a={c[0]:.3f}$, $\\delta={c[1]}$" for c in pj.columns]
    lam_vals = sorted(set(c[0] for c in pj.columns)); del_vals = sorted(set(c[1] for c in pj.columns))
    body = ("$\\Delta$ & $\\eta$ & " + " & ".join(f"\\multicolumn{{3}}{{c}}{{Journal, $\\lambda_a={l:.3f}$}}" for l in lam_vals) + " & \\multicolumn{3}{c}{Conference} \\\\\n"
            + " & & " + " & ".join(" & ".join(f"$\\delta$={d}" for d in del_vals) for _ in lam_vals) + " & " + " & ".join(f"$\\delta$={d}" for d in pc.columns) + " \\\\\n\\midrule\n" + "\n".join(rows))
    note = "Journal scenarios with $L=0.3$. $T^*$ is searched on 3--60 days, so 60 denotes a corner at or beyond 60 days."
    (OUT / "tabS5_sensitivity.tex").write_text(wrap(body, "Sensitivity of $T^*$ to the quality penalty $\\Delta$, the compression elasticity $\\eta$, the acceptance speed $\\lambda_a$ and the delay cost $\\delta$", "tab:sensitivity", "llrrrrrrrrrrrr", note, "\\scriptsize", colsep="2.5pt", env="sidewaystable", width="0.95\\textheight"))


if __name__ == "__main__":
    table1(); table2(); table3(); table4(); table5(); table6(); or2_tables()
    print("tablas escritas en", OUT, sorted(p.name for p in OUT.glob("*.tex")))
