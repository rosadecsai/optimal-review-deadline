"""Builds the LaTeX tables of the revised manuscript and of Online Resource 3
from output/results.json (original pipeline) and revision_results.json."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

REPO, REV, TAB = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
RES = json.loads((REPO / "output" / "results.json").read_text())
R = json.loads((REV / "revision_results.json").read_text())
FC = json.loads((REV / "final_checks.json").read_text())      # checks of the final version (code/final_checks.py)
TAB.mkdir(parents=True, exist_ok=True)


def f(x, d=3):
    s = f"{x:.{d}f}"
    return s.replace("-", "$-$") if s.startswith("-") else s


def ci(lo, hi, d=3):
    return f"[{f(lo, d)}, {f(hi, d)}]"


# ------------------------------------------------ Table 5: Parkinson slopes ---
a1 = {r["stat"]: r for r in R["A1"]["table"]}
boot = {r["stat"]: r for r in RES["E3"]["slopes"]}
rows = [("s_q25", "$Q_{25}$, days before deadline", "0"), ("s_q50", "Median, days before deadline", "0"),
        ("s_q75", "$Q_{75}$, days before deadline", "0"), ("s_q90", "$Q_{90}$, days before deadline", "0"),
        ("x_q50", "Median, days since assignment", "1"), ("x_q50_all", "\\rev{\\quad same, all reviews}", "\\rev{1}"),
        ("share_atom", "Share in final 48 h \\rev{(points)}", "--"),
        ("share_late", "Share late \\rev{(points)}", "--"), ("gamma", "Anchoring rate $\\hat\\gamma$", "--")]
body = []
for key, lab, pk in rows:
    if key == "x_q50":
        s = a1["s_q50"]; b = boot["x_q50"]
        body.append(f"{lab} & {f(1 - s['slope'])} & {ci(1 - s['hi'], 1 - s['lo'])} & \\rev{{{s['p_perm']:.3f}$^{{\\dagger}}$}} & {ci(b['lo'], b['hi'])} & {pk} \\\\")
        continue
    if key == "x_q50_all":
        s = FC["F3"]["all_reviews"]["q50"]
        body.append(f"{lab} & \\rev{{{f(s['passthrough'])}}} & \\rev{{{ci(s['lo'], s['hi'])}}} & \\rev{{{s['p_perm']:.3f}$^{{\\dagger}}$}} & -- & {pk} \\\\")
        continue
    s = a1[key]; b = boot[key]; k100 = 100 if key in ("share_atom", "share_late") else 1; d = 2 if k100 == 100 else 3
    body.append(f"{lab} & {f(k100 * s['slope'], d)} & {ci(k100 * s['lo'], k100 * s['hi'], d)} & {s['p_perm']:.3f} & {ci(k100 * b['lo'], k100 * b['hi'], d)} & {pk} \\\\")
(TAB / "tab3_parkinson.tex").write_text(r"""\begin{table}[t]
\caption{Delivery quantiles and the review window: slopes across the seven cohorts with a documented window}\label{tab:parkinson}
\footnotesize
\setlength{\tabcolsep}{2pt}
\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lrlrlr@{}}
\toprule
Statistic & Slope & 95\% CI (cohorts) & Perm.\ $p$ & Within-cohort & \rev{Full anchoring} \\
\midrule
""" + "\n".join(body) + r"""
\bottomrule
\end{tabular*}
\footnotetext{\rev{Slope from an ordinary least-squares (OLS) regression} of the cohort statistic on the effective window $T$ across the seven cohorts 2020--2026 ($T=12$--21 days)\rev{; the slopes of the two shares are in percentage points per day of window. All statistics are over on-time reviews, except the shares and the row marked ``all reviews'', which include late ones}. 95\% CI (cohorts): the seven cohorts are the units of observation (Student $t$, 5 degrees of freedom). Perm.\ $p$: exact two-sided permutation test of a zero slope over the $7!$ assignments of the cohort statistics to the windows. \rev{$^{\dagger}$For the two rows in days since assignment the test is of a slope of one, that is, of a zero slope of days before the deadline.} Within-cohort: interval from manuscripts resampled within each cohort (500 replications); it measures only the sampling precision of each cohort statistic, not the uncertainty of a slope estimated from seven cohorts, and is reported for comparison with the submitted version. \rev{Full anchoring: the limit of Proposition~\ref{prop:invariance} in which nothing is delivered in the days that a longer window adds (days before the deadline do not depend on $T$; days since assignment move one-for-one); the model itself implies small positive slopes of days before the deadline (Online Resource~3, Table~S27).} The 12-day cohort (2023) has leverage 0.77\rev{, the weight of a single observation in a fitted slope; values close to one mean that the slope is largely decided by that observation}.}
\end{table}
""")

# ----------------------------------------------------- Table 9: T* panels ----
ts = pd.DataFrame(RES["E7"]["tstar"])


def row(r, lead):
    mass = 100 * (r["share_atom_Tstar"] + r["share_late_Tstar"])
    losses = " & ".join(f"{r[f'loss_T{T}']:.1f}" for T in (7, 14, 21, 28, 42))
    return f"{lead} & {r['delta']:g} & {r['eta']:.1f} & {int(r['T_star'])} & {r['a_Tstar']:.2f} & {mass:.0f} & {losses} \\\\"


pa = [row(r, "--") for _, r in ts[ts["conference"]].sort_values(["eta", "delta"]).iterrows()]
pb = [row(r, f"{r['L']:.1f}") for _, r in ts[~ts["conference"]].sort_values(["eta", "L", "delta"]).iterrows()]
(TAB / "tab6_tstar.tex").write_text(r"""\begin{sidewaystable}[p]
\caption{The window that maximises the editor's expected value $V(T)$, by scenario. Panel A uses the ICLR delivery profile; Panel B is a model-based illustration for journals built on parameters taken from outside our data}\label{tab:tstar}
\scriptsize
\setlength{\tabcolsep}{3pt}
\begin{tabular*}{0.95\textheight}{@{\extracolsep{\fill}}lllrrrrrrrr@{}}
\toprule
$L$ & $\delta$ & $\eta$ & $T^*$ (days) & $a(T^*)$ & Deadline mass (\%) & \multicolumn{5}{c}{Loss at $T$ = 7 / 14 / 21 / 28 / 42 days (\% of $\bar q$)} \\
\midrule
\multicolumn{11}{@{}p{0.93\textheight}}{\textbf{Panel A. Conference}: reviewers committed \rev{in advance, for the conference as a whole} ($a\equiv1$, $L=0$)} \\
""" + "\n".join(pa) + r"""
\midrule
\multicolumn{11}{@{}p{0.93\textheight}}{\textbf{Panel B. Journal (illustrative)}: acceptance curve $a(T)=\bar a(1-e^{-\lambda_a T})$ through the two arms of \citet{chetty2014}, $\bar a=0.684$, $\lambda_a=0.099$} \\
""" + "\n".join(pb) + r"""
\bottomrule
\end{tabular*}
\footnotetext{How to read the table. Each row is one combination of the parameters that we set (Table~\ref{tab:inputs}): $L$, the loss when an invitation is declined; $\delta$, the cost of one day of delay; $\eta$, how fast the last-minute shortfall grows as the window shortens, $\Delta(T)=\Delta\,(20.5/T)^{\eta}$ \rev{with $\Delta=0.085$ for reports delivered in the final 48 hours and 0.095 for late reports, and 20.5 days the mean window of the cohorts on which they are measured} ($\eta=0$: the shortfall does not depend on the window). All three are in units of the value $\bar q=1$ of a complete report. $T^*$: the window with the highest $V(T)$ on 3--60 days \rev{(3--90 in Sect.~\ref{sec:joint} and Online Resource~3)}. $a(T^*)$: probability that the invitation is accepted at $T^*$. Deadline mass: predicted share of reviews delivered in the final 48 hours or late at $T^*$. Loss at $T$: $100\,[V(T^*)-V(T)]$, the value forgone by choosing $T$ instead of $T^*$. Delivery moments for every $T$ come from the 2024--2026 hazard profile (Sect.~\ref{sec:calib})\rev{, estimated on windows of 20--21 days: the columns for 7 days and for 28 and 42 days, and every optimum below 12 or above 21 days, are extrapolations}. Panel B does not use journal data: it combines the conference delivery profile with an acceptance curve fitted to \rev{the experiment of \citet{chetty2014}} at one economics journal, and its numbers move with the assumptions examined in Sect.~\ref{sec:joint}.}
\end{sidewaystable}
""")

# --------------------------------------------------------- Online Resource 3 -
P = []


def table(caption, label, header, lines, note, colspec):
    P.append("\\begin{table}[htbp]\n\\caption{" + caption + "}\\label{" + label + "}\n\\centering\\footnotesize\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{" + colspec + "}\n\\toprule\n"
             + header + " \\\\\n\\midrule\n" + "\n".join(lines) + "\n\\bottomrule\n\\end{tabular}\n\n\\smallskip\n\\begin{minipage}{0.95\\textwidth}\\footnotesize " + note + "\\end{minipage}\n\\end{table}\n")


# S6 cohort statistics and slopes
e3 = pd.DataFrame(RES["E3"]["table"]).sort_values("year")
lines = [f"{int(r.year)} & {int(r['T'])} & {r['s_q25']:.2f} & {r['s_q50']:.2f} & {r['s_q75']:.2f} & {r['s_q90']:.2f} & {100 * r['share_atom']:.1f} & {100 * r['share_late']:.1f} & {r['gamma']:.3f} \\\\" for _, r in e3.iterrows()]
lines.append("\\midrule")
for lab, k in (("Slope on $T$", "slope"), ("Standard error", "se"), ("95\\% CI, lower", "lo"), ("95\\% CI, upper", "hi"), ("Permutation $p$", "p_perm"), ("Leave-one-out, min", "loo_min"), ("Leave-one-out, max", "loo_max")):
    v = [a1[s][k] for s in ("s_q25", "s_q50", "s_q75", "s_q90", "share_atom", "share_late", "gamma")]
    v = [x * (100 if s in (4, 5) and k != "p_perm" else 1) for s, x in enumerate(v)]
    lines.append(f"{lab} & & " + " & ".join(f(x, 3) for x in v) + " \\\\")
table("Cohort statistics of on-time delivery and cohort-level inference on their slopes with respect to the window", "tab:s6",
      "Cohort & $T$ & $Q_{25}$ & Median & $Q_{75}$ & $Q_{90}$ & Final 48 h (\\%) & Late (\\%) & $\\hat\\gamma$", lines,
      "Quantiles are in days before the deadline, over on-time reviews. Slopes are OLS across the seven cohorts; the interval uses Student's $t$ with 5 degrees of freedom; the permutation $p$-value is exact (5,040 orderings); leave-one-out gives the smallest and largest slope when one cohort is dropped. Shares and their slopes are in percentage points. The correlation between $T$ and the calendar year is "
      + f"{R['A1']['corr_T_year']:.2f}.", "lrrrrrrrr")

# S7 acceptance forms
a3 = pd.DataFrame(R["A3"]["table"])
names = [("exponential (submitted)", "Saturating exponential (baseline)"), ("hyperbolic", "Hyperbolic, $\\bar a\\,T/(T+k)$"), ("logarithmic", "Logarithmic, $c+b\\ln T$"), ("linear (chord)", "Linear (chord through the two arms)")]
lines = []
for key, lab in names:
    x = a3[a3["form"] == key]
    g = lambda d, e, c: x[(x["delta"] == d) & (x["eta"] == e)][c].iloc[0]
    lines.append(f"{lab} & {g(0.005, 0.0, 'a7'):.2f} & {g(0.005, 0.0, 'a14'):.2f} & {g(0.005, 0.0, 'a21'):.2f} & " + " & ".join(str(int(g(d, e, 'T_star'))) for e in (0.0, 1.3) for d in (0.0025, 0.005, 0.01)) + " \\\\")
table("Journal illustration under four acceptance curves that pass through the same two experimental arms, $a(28)=0.641$ and $a(45)=0.676$", "tab:s7",
      "Acceptance curve & $a(7)$ & $a(14)$ & $a(21)$ & \\multicolumn{3}{c}{$T^*$, $\\eta=0$} & \\multicolumn{3}{c}{$T^*$, $\\eta=1.3$} \\\\\n & & & & $\\delta$=0.0025 & 0.005 & 0.01 & 0.0025 & 0.005 & 0.01", lines,
      "$L=0.3$, $\\Delta=0.085$, ICLR delivery profile, $T^*$ searched on 3--90 days (90 denotes a corner at or beyond 90 days). Every increasing concave acceptance curve through the two arms lies on or below the chord outside the interval 28--45 days, so the chord is the curve under which shortening the window costs least acceptance.", "lrrrrrrrrr")

# S8 curvature with anchored level; pass-through
lam = pd.DataFrame(R["A3_lambda"])
lines = [f"{l:.3f} & {lam[lam['lam_a'] == l]['abar'].iloc[0]:.3f} & {lam[lam['lam_a'] == l]['a14'].iloc[0]:.2f} & " + " & ".join(str(int(lam[(lam['lam_a'] == l) & (lam['delta'] == d)]['T_star'].iloc[0])) for d in (0.0025, 0.005, 0.01)) + " \\\\" for l in sorted(lam["lam_a"].unique())]
table("Journal illustration by the curvature $\\lambda_a$ of the acceptance curve, with its level anchored at the 28-day arm, $a(28)=0.641$", "tab:s8",
      "$\\lambda_a$ & $\\bar a$ & $a(14)$ & $T^*$, $\\delta$=0.0025 & 0.005 & 0.01", lines,
      "$L=0.3$, $\\Delta=0.085$, $\\eta=0$. In the submitted Table~S5 the asymptote $\\bar a$ was held at 0.684 while $\\lambda_a$ varied, which also moved $a(28)$ away from the experimental value; here the level is kept at the arm closest to the optimum.", "rrrrrr")

pc = pd.DataFrame(R["A4"]["closed"])
lines = []
for rho in sorted(pc["rho"].unique(), reverse=True):
    for L in (0.1, 0.3, 0.6):
        x = pc[(pc["rho"] == rho) & (pc["L"] == L)]
        lines.append(f"{rho:.2f} & {L:.1f} & " + " & ".join(f"{x[x['delta'] == d]['T_closed'].iloc[0]:.1f}" for d in (0.0025, 0.005, 0.01)) + " \\\\")
table("Journal illustration with delivery anchored on the experiment of Chetty et al.\\ (2014) instead of the conference profile: closed form with partial pass-through $\\rho$", "tab:s9",
      "$\\rho$ & $L$ & $T^*$, $\\delta$=0.0025 & 0.005 & 0.01", lines,
      "Expected delivery $\\mathbb{E}[\\tau\\mid T]=36+\\rho\\,(T-28)$ days (the median of the four-week arm), no shortfall ($\\Delta=0$), baseline acceptance curve. $\\rho=0.71$ is the ratio of the 12-day fall in median turnaround to the 17-day cut in the deadline in that experiment; $\\rho=0.93$ is the slope of the median delivery date on the window across the ICLR cohorts; $\\rho=1$ is full anchoring. Closed form of Proposition~4 as generalised in Remark~2 of the main text; it agrees with a numerical maximisation to within "
      + f"{R['A4']['max_diff']:.4f} days.", "rrrrr")

# S10 tails
a5 = pd.DataFrame(R["A5"]["table"]).set_index("tail")
tn = {"power": "Power law (baseline)", "exponential": "Exponential", "truncated": "Power law, cut at 30 d", "floor": "Constant beyond 15 d"}
lines = [f"{tn[t]} & {100 * a5.loc[t, 'reach48_T21']:.0f} & {100 * a5.loc[t, 'reach48_T45']:.0f} & {100 * a5.loc[t, 'reach48_T400']:.0f} & "
         + " & ".join(str(int(a5.loc[t, c])) for c in ("T_conf_eta1.3_0.0025", "T_conf_eta1.3_0.005", "T_conf_eta1.3_0.01", "T_jour_eta0_0.0025", "T_jour_eta0_0.005", "T_jour_eta0_0.01")) + " \\\\" for t in tn]
table("Sensitivity to the extrapolation of the hazard beyond the observed range (more than 15 days before the deadline)", "tab:s10",
      "Tail of the hazard & \\multicolumn{3}{c}{Unwritten at 48 h (\\%)} & \\multicolumn{3}{c}{$T^*$ conference, $\\eta=1.3$} & \\multicolumn{3}{c}{$T^*$ journal, $\\eta=0$} \\\\\n & $T$=21 & 45 & 400 & $\\delta$=0.0025 & 0.005 & 0.01 & 0.0025 & 0.005 & 0.01", lines,
      f"``Unwritten at 48 h'': share of reviewers who reach the final 48 hours with the report undelivered. Baseline: $h_s={R['A5']['pl_a']:.2f}\\,(1+s)^{{-{R['A5']['pl_b']:.2f}}}$ for $s\\ge15$. Exponential: $h_s=\\exp({R['A5']['exp_c']:.2f}-{R['A5']['exp_k']:.3f}\\,s)$, fitted on the same days ($s=2$--15). "
      "The third variant sets the hazard to zero more than 30 days before the deadline; the fourth holds it at its value 15 days before (the hazard never stops, so with a long enough window almost every report is written early). Journal: $L=0.3$, baseline acceptance curve.", "lrrrrrrrrr")

# S11 Delta
a6 = pd.DataFrame(R["A6"])
lines = []
for D in sorted(a6["Delta"].unique()):
    for eta in (0.0, 1.3):
        x = a6[(a6["Delta"] == D) & (a6["eta"] == eta)]
        lines.append(f"{D:.3f} & {eta:.1f} & " + " & ".join(str(int(x[x['delta'] == d][c].iloc[0])) for c in ("T_conference", "T_journal") for d in (0.0025, 0.005, 0.01)) + " \\\\")
table("The window that maximises $V(T)$ as the shortfall $\\Delta$ assigned to last-minute and late reports goes from zero to 0.30", "tab:s11",
      "$\\Delta$ & $\\eta$ & \\multicolumn{3}{c}{$T^*$ conference} & \\multicolumn{3}{c}{$T^*$ journal ($L=0.3$)} \\\\\n & & $\\delta$=0.0025 & 0.005 & 0.01 & 0.0025 & 0.005 & 0.01", lines,
      "$\\Delta=0$ corresponds to the reading that a shorter report is not a less valuable one. $T^*$ searched on 3--90 days; 3 denotes a corner at the shortest window considered. One value of $\\Delta$ is applied here to reports of the final 48 hours and to late reports alike; the scenario table of the manuscript uses 0.085 for the former and 0.095 for the latter, which is why its conference optima at $\\eta=1.3$ are 25, 19 and 14 days where the row $\\Delta=0.085$ of this table gives 24, 18 and 14. With $\\eta=0$ the journal optimum falls by a day as $\\Delta$ rises: a larger shortfall lowers the value of every accepted report, and with it the marginal value of acceptance, by more than it raises the gain from fewer last-minute reports (Proposition~5 of the manuscript).", "rrrrrrrr")

# S12 joint
a7 = R["A7"]
qs = ["5", "25", "50", "75", "95"]
lines = [f"Conference, all draws & " + " & ".join(f"{a7['conference']['q'][q]:.0f}" for q in qs) + " \\\\",
         f"Conference, $\\eta<0.5$ & " + " & ".join(f"{a7['conference']['eta<0.5'][q]:.0f}" for q in qs) + " \\\\",
         f"Conference, $\\eta\\ge1$ & " + " & ".join(f"{a7['conference']['eta>=1'][q]:.0f}" for q in qs) + " \\\\",
         f"Journal, all draws & " + " & ".join(f"{a7['journal']['q'][q]:.0f}" for q in qs) + " \\\\"]
mp = a7["journal_map"]
for ln, ll in (("low", "low"), ("mid", "middle"), ("high", "high")):
    for dn, dl in (("low", "low"), ("mid", "middle"), ("high", "high")):
        v = mp[f"lam_{ln}|delta_{dn}"]
        lines.append(f"Journal, $\\lambda_a$ {ll}, $\\delta$ {dl} & {v[0]:.0f} & & {v[1]:.0f} & & {v[2]:.0f} \\\\")
table("Joint sensitivity: percentiles of the window that maximises $V(T)$ over " + f"{a7['n']:,}".replace(",", "{,}") + " Latin-hypercube draws of all parameters at once", "tab:s12",
      "Scenario & 5th & 25th & Median & 75th & 95th", lines,
      "Ranges, sampled independently: $\\delta$ log-uniform on 0.0025--0.01; $L$ uniform on 0.1--0.6; $\\Delta$ uniform on 0--0.30; $\\eta$ uniform on 0--2; $\\lambda_a$ log-uniform on 0.05--0.20 with $a(28)=0.641$; pass-through $\\rho$ uniform on 0.7--1 (journal); late hazard $\\mu$ uniform on 0.20--0.42; tail of the hazard drawn from the four variants of Table~\\ref{tab:s10}. "
      f"The ranges are not beliefs: the percentiles describe how the output moves over the box, not a probability. Terciles of $\\lambda_a$: below {a7['lam_terciles'][0]:.3f}, {a7['lam_terciles'][0]:.3f}--{a7['lam_terciles'][1]:.3f}, above; of $\\delta$: below {a7['delta_terciles'][0]:.4f}, {a7['delta_terciles'][0]:.4f}--{a7['delta_terciles'][1]:.4f}, above. "
      f"$V(T)$ is single-peaked on 3--90 days in {100 * a7['conference']['unimodal']:.0f}\\% of the conference draws and {100 * a7['journal']['unimodal']:.0f}\\% of the journal draws. The 25th and 75th percentiles were not computed for the rows by tercile. The journal rows use the exponential family of acceptance curves only; the functional form is varied in Table~\\ref{{tab:s7}} and in the last section of this resource.", "lrrrrr")

nm = {"lam_a": "Acceptance curvature $\\lambda_a$", "delta": "Delay cost $\\delta$", "eta": "Compression elasticity $\\eta$", "Delta": "Shortfall $\\Delta$", "L": "Loss from a decline $L$",
      "rho": "Pass-through $\\rho$", "tail": "Tail of the hazard", "mu": "Late hazard $\\mu$", "interactions": "Interactions (remainder)"}
lines = [f"{nm[k]} & " + (f"{100 * a7['S_conference'][k]:.1f}" if k in a7["S_conference"] else "--") + f" & {100 * a7['S_journal'][k]:.1f} \\\\" for k in nm]
table("First-order sensitivity indices: share of the variance of $T^*$ across the joint draws explained by each parameter on its own (\\%)", "tab:s13",
      "Parameter & Conference & Journal", lines,
      "Var$(\\mathbb{E}[T^*\\mid x])/$Var$(T^*)$, estimated with 20 equal-count bins of each continuous parameter and the four categories of the tail; the remainder is the share due to interactions.", "lrr")

# S14 before/after
a9 = R["A9"]
on = {"log_words": "log words", "confidence": "Confidence", "score": "Overall score", "dev_panel": "Deviation from co-reviewers"}
lines = [f"{on[k]} & {f(a9[k]['pre']['est'])} ({a9[k]['pre']['se']:.3f}) & {f(a9[k]['post']['est'])} ({a9[k]['post']['se']:.3f}) & {f(a9[k]['post_excl_2025']['est'])} ({a9[k]['post_excl_2025']['se']:.3f}) & {f(a9[k]['diff'])} & {a9[k]['z']:.2f} \\\\" for k in on]
table("Final-24-hour coefficients before and after 2024: inverse-variance pooling of the cohort-by-cohort within-manuscript estimates", "tab:s14",
      "Outcome & 2018--2023 & 2024--2026 & 2024 and 2026 & Difference & $z$", lines,
      "Pooled from the cohort estimates of Table~S3 (reference: more than four days before the deadline); standard errors in parentheses. ICLR 2024 is the first cohort reviewed after the public release of ChatGPT and also the first with the four-field review form; the two cannot be separated. The 2025 cohort hosted a randomised experiment of LLM feedback to reviewers, so the fourth column drops it. "
      f"Without the 12-day cohort of 2023 the pooled length coefficient of the earlier cohorts is {f(R['A12']['log_words_pre_excl_2023']['est'])} ({R['A12']['log_words_pre_excl_2023']['se']:.3f}), still larger than that of 2024--2026 ($z={R['A12']['log_words_pre_excl_2023']['z_vs_post']:.2f}$).", "lrrrrr")

# S15 kappa (manuscript-clustered bootstrap from the review-level file; falls back to the approximation of A10)
_kb = REV / "kappa_bootstrap.json"
if _kb.exists():
    KB = json.loads(_kb.read_text()); k = pd.DataFrame(KB["table"]); pl = KB["pooled"]
    dn = "early-early minus late-late"
    lines = [f"{int(r.year)} & {r['early-early']:.2f} ({r['se_early-early']:.2f}) & {int(r['n_early-early']):,} & {r['atom-atom']:.2f} ({r['se_atom-atom']:.2f}) & {r['late-late']:.2f} ({r['se_late-late']:.2f}) & {int(r['n_late-late']):,} & {f(r[dn], 2)} & {f(r['lo_' + dn], 2)} to {f(r['hi_' + dn], 2)} \\\\".replace(",", "{,}") for _, r in k.iterrows()]
    x1, x2, x3 = pl["early-early minus late-late | 2021-2025"], pl["early-early minus atom-atom | 2021-2025"], pl["atom-atom minus late-late | 2021-2025"]
    y1 = pl["early-early minus late-late | 2018-2025"]
    table("Cohen's $\\kappa$ on the accept vote by the timing of the two reviews, with standard errors from a bootstrap over manuscripts", "tab:s15",
          "Cohort & Early--early & Pairs & Final 48 h--final 48 h & Late--late & Pairs & Early $-$ late & 95\\% interval", lines,
          f"Standard errors in parentheses and intervals from {KB['n_boot']:,} bootstrap resamples of manuscripts within each cohort (every $\\kappa$ is a function of sums over manuscripts, so pairs that share a manuscript are resampled together).".replace(",", "{,}")
          + f" Pooled over 2021--2025 (the cohorts with at least 225 early--early pairs) with inverse-variance weights, $\\kappa$ of early--early pairs exceeds that of late--late pairs by {x1['est']:.3f} (SE {x1['se']:.3f}; 95\\% interval {x1['lo']:.3f} to {x1['hi']:.3f}; positive in {x1['n_positive']} of {x1['k']} cohorts) and that of final-48-hour pairs by {x2['est']:.3f} ({x2['se']:.3f}; {x2['lo']:.3f} to {x2['hi']:.3f}; positive in {x2['n_positive']} of {x2['k']}); final-48-hour pairs exceed late--late pairs by {x3['est']:.3f} ({x3['se']:.3f}; {x3['lo']:.3f} to {x3['hi']:.3f})."
          + f" The unweighted means of the cohort differences are {x1['unweighted']:.3f} ({x1['unweighted_lo']:.3f} to {x1['unweighted_hi']:.3f}) and {x2['unweighted']:.3f} ({f(x2['unweighted_lo'], 3)} to {x2['unweighted_hi']:.3f}). Over all eight cohorts the first pooled difference is {y1['est']:.3f} ({y1['lo']:.3f} to {y1['hi']:.3f}); the cohorts of 2018--2020 have 29 to 106 early--early pairs. The archived values of $\\kappa$ are reproduced exactly from the review-level file (\\texttt{{code/kappa\\_bootstrap.py}}). Early: four days or more before the deadline; pairs that involve a review delivered two to four days before the deadline (four of the ten types of pair) are not shown.", "lrrrrrrr")
else:
    k = pd.DataFrame(R["A10"]["table"])
    lines = [f"{int(r.year)} & {r['early-early']:.2f} ({r['se_early-early']:.2f}) & {int(r['n_early-early']):,} & {r['atom-atom']:.2f} ({r['se_atom-atom']:.2f}) & {r['late-late']:.2f} ({r['se_late-late']:.2f}) & {int(r['n_late-late']):,} & {f(r['early-early'] - r['late-late'], 2)} \\\\".replace(",", "{,}") for _, r in k.iterrows()]
    pl = R["A10"]["pooled"]
    x1, x2 = pl["early-early minus late-late | 2021-2025"], pl["early-early minus atom-atom | 2021-2025"]
    table("Cohen's $\\kappa$ on the accept vote by the timing of the two reviews, with approximate standard errors", "tab:s15",
          "Cohort & Early--early & Pairs & Final 48 h--final 48 h & Late--late & Pairs & Early $-$ late", lines,
          "Approximate standard errors in parentheses, treating pairs as independent and using the cohort's accept rate for the marginals; pairs that share a manuscript are not independent, so these errors are lower bounds. "
          f"Pooled over 2021--2025 (the cohorts with at least 225 early--early pairs), $\\kappa$ of early--early pairs exceeds that of late--late pairs by {x1['est']:.3f} (approximate SE {x1['se']:.3f}; positive in {x1['n_positive']} of {x1['k']} cohorts) and that of final-48-hour pairs by {x2['est']:.3f} ({x2['se']:.3f}; positive in {x2['n_positive']} of {x2['k']}).", "lrrrrrr")

# S16 reviews per manuscript
a11 = pd.DataFrame(R["A11"])
lines = [f"{int(r.year)} & {int(r.manuscripts):,} & {r.mean_reviews:.2f} & {int(r.modal)} & {100 * r.share_lt3:.1f} & {100 * r.share_ge4:.1f} \\\\".replace(",", "{,}") for r in a11.itertuples()]
table("Official reviews per reviewed manuscript, by cohort", "tab:s16", "Cohort & Manuscripts & Mean reviews & Mode & Fewer than three (\\%) & Four or more (\\%)", lines,
      "Manuscripts with at least one official review, from the paper-level file of the repository; the counts differ from those of the cohort table of the main text by at most 12 manuscripts (2018), because a few reviews are not linked to a paper record or fall outside the analysis sample. Assignments, declines and withdrawals of reviewers are not exposed by the API, so non-completion is not observed directly; the table shows that almost every reviewed manuscript ends with at least three reviews, including those supplied by replacement reviewers after the deadline.", "rrrrrr")

head = r"""\documentclass[11pt]{article}
\usepackage[a4paper,margin=2.2cm]{geometry}
\usepackage{booktabs,amsmath,amssymb,rotating,url}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\renewcommand{\thetable}{S\arabic{table}}
\setcounter{table}{7}
\begin{document}
\begin{center}
{\large\bfseries Online Resource 3\\[0.5ex] Analyses added in the revision: sensitivity of the calibration, two review windows inside one journal, and final checks}\\[1.2ex]
Supplement to: \emph{How Long Should Reviewers Get? Parkinson's Law and the Optimal Review Deadline in 209,000 Timestamped Conference Reviews}\\[0.6ex]
\emph{Scientometrics}\\[0.4ex]
Jos\'e A. Garc\'ia, Rosa Rodriguez-S\'anchez, J. Fdez-Valdivia\\
Universidad de Granada, Spain. Corresponding author: jags@decsai.ugr.es
\end{center}

\noindent This resource reports the analyses added in the revision of the manuscript. Its tables are numbered after those of Online Resource~1 (Tables~S1--S5) and Online Resource~2 (Tables~S6--S7). The first section covers the conference data: cohort-level inference on the effect of the window (Table~\ref{tab:s6}), the journal illustration under alternative acceptance curves, curvatures and delivery assumptions (Tables~\ref{tab:s7}--\ref{tab:s9}), alternative tails of the hazard (Table~\ref{tab:s10}), the shortfall $\Delta$ from zero (Table~\ref{tab:s11}), the joint sensitivity analysis (Tables~\ref{tab:s12}--\ref{tab:s13}), timing coefficients before and after 2024 (Table~\ref{tab:s14}), agreement by timing pair with standard errors from a bootstrap over manuscripts (Table~\ref{tab:s15}) and the number of reviews per manuscript (Table~\ref{tab:s16}). Every number of these tables is produced by \texttt{code/revision\_analyses.py} from the archived outputs of the original pipeline (the bootstrap of the agreement table, by \texttt{code/kappa\_bootstrap.py} from the review-level file); the script first reproduces the 24 optima of the main table exactly. The second section reports the within-journal comparison of two review windows at TMLR (Tables~\ref{tab:tmlr-acc}--\ref{tab:tmlr-sel}), and the third the checks made on the final version (Tables~\ref{tab:fc-haz}--\ref{tab:fc-opt}).

\section*{Conference data: cohort-level inference and sensitivity of the calibration}

"""
# ----------------------------------------------- checks of the final version (code/final_checks.py)
P_main = list(P); P.clear()
F1, F2, F3, F4, F5, F6, F7, F8 = (FC[k] for k in ("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8"))
pcd = lambda v, d=1: f"{100 * v:.{d}f}"
lines = []
for d_ in range(0, 12):
    r = F1["per_day"][str(d_)]; lo_, hi_ = F1["hazard_range_by_day"][str(d_)]
    lines.append(f"{d_} & {pcd(lo_)} & {pcd(hi_)} & {f(r['slope'])} & {f(r['lo'])} to {f(r['hi'])} \\\\")
lines.append("\\midrule")
for lab, r in (("Days 0--3, pooled", F1["near"]), ("Days 4--11, pooled", F1["far"]), ("Days 0--11, pooled", F1["pooled_log_hazard_slope"])):
    lines.append(f"{lab} & & & {f(r['slope'])} & {f(r['lo'])} to {f(r['hi'])} \\\\")
table("ICLR: the daily hazard on the days before the deadline that all seven cohorts share, and its slope on the window", "tab:fc-haz",
      "Days before the deadline & \\multicolumn{2}{c}{Hazard across cohorts (\\%)} & Slope of the log hazard & 95\\% interval \\\\\n & lowest & highest & on the window $T$ &", lines,
      "Proposition~1 of the manuscript implies that the hazard on a given day before the deadline does not depend on the window, that is, a zero slope in every row. Slope: ordinary least squares, across the seven cohorts with a documented window (12--21 days), of the logarithm of the daily hazard on the window; intervals from Student's $t$ with 5 degrees of freedom. Pooled rows: the log hazard of each cohort relative to the mean of the seven cohorts on that day, averaged over the days. "
      f"The pooled slope of {f(F1['pooled_log_hazard_slope']['slope'])} corresponds to {F1['pooled_pct_per_day']:.1f}\\% more hazard per day of window, with an interval from {f(F1['pooled_pct_lo'], 1)}\\% to {F1['pooled_pct_hi']:.1f}\\%; {F1['n_days_slope_significant']} of the 12 daily slopes differs from zero at the 5\\% level. The test has little power, and the cohorts with the longer windows are also the recent, more punctual ones. Hazards from \\texttt{{output/tables/hazard\\_by\\_day.csv}}.", "rrrrr")
qn = {"q25": "$Q_{25}$", "q50": "Median", "q75": "$Q_{75}$", "q90": "$Q_{90}$", "share_final48h": "Final 48 h (points)", "share_late": "Late (points)"}
a1k = {"q25": "s_q25", "q50": "s_q50", "q75": "s_q75", "q90": "s_q90", "share_final48h": "share_atom", "share_late": "share_late"}
lines = []
for k, lab in qn.items():
    k100 = 100 if k.startswith("share") else 1; dd = 2 if k100 == 100 else 3
    e_ = a1[a1k[k]]; pb_ = F5["published"][k]
    imp = " & ".join(f(F2[t_][k]["implied_slope"], dd) for t_ in ("power", "exponential", "floor"))
    lines.append(f"{lab} & {f(k100 * e_['slope'], dd)} & {f(k100 * e_['lo'], dd)} to {f(k100 * e_['hi'], dd)} & {f(pb_['slope'], dd)} & {f(pb_['lo'], dd)} to {f(pb_['hi'], dd)} & {imp} \\\\")
lines.append("\\midrule")
lines.append("\\multicolumn{8}{l}{\\emph{Pass-through of the window to the median delivery date}} \\\\")
pt_ = F3["on_time"]; pa_ = F3["all_reviews"]["q50"]; p6_ = F3["on_time_without_2023"]; pp_ = F5["published"]["passthrough_median"]; a6_ = F3["all_reviews"]["q50"]["without_2023"]
lines.append(f"On-time reviews & {f(pt_['passthrough'])} & {f(pt_['lo'])} to {f(pt_['hi'])} & {f(pp_['est'])} & {f(pp_['lo'])} to {f(pp_['hi'])} & & & \\\\")
lines.append(f"All reviews & {f(pa_['passthrough'])} & {f(pa_['lo'])} to {f(pa_['hi'])} & & & & & \\\\")
lines.append(f"On-time, without 2023 & {f(p6_['passthrough'])} & {f(p6_['lo'])} to {f(p6_['hi'])} & & & & & \\\\")
lines.append(f"All, without 2023 & {f(a6_['passthrough'])} & {f(a6_['lo'])} to {f(a6_['hi'])} & & & & & \\\\")
table("ICLR: the slopes of the cohort table of the manuscript under the published deadlines and windows, and the slopes that the estimated hazards imply", "tab:fc-slopes",
      "Statistic & \\multicolumn{2}{c}{Effective deadlines} & \\multicolumn{2}{c}{Published deadlines} & \\multicolumn{3}{c}{Implied by the profile} \\\\\n & slope & 95\\% interval & slope & 95\\% interval & power & expon. & const.",
      lines,
      "Slopes across the seven cohorts with a documented window. Quantiles: days before the deadline, over on-time reviews; shares: over all reviews, in percentage points per day of window; pass-through: days since assignment per day of window, with and without late reviews and with and without the 12-day cohort. Effective deadlines: as in the manuscript (windows of 20, 16, 19, 12, 21, 20 and 21 days). Published deadlines: the deadlines and windows of the reviewer guides (20, 16, 18, 12, 21, 20 and 22 days), which differ from the effective ones by one day in 2022 and in 2026; in 2022 the reviews of the day we take as an extension then count as late. "
      "Implied: the slopes obtained when the same hazard profile (2024--2026) is applied to the seven windows, as Proposition~1 prescribes, under three tails of the hazard beyond 15 days. The model does not imply zero slopes of the quantiles: a longer window adds early days on which some reports are delivered. The observed slopes are two to three times the implied ones at the lower quartile and the median and close to them at the ninth decile; none of the differences is distinguishable from zero with seven cohorts. "
      f"The last rows give the pass-through of the median; the lower end of its interval is {pt_['lo']:.2f} with the effective deadlines and on-time reviews, {pa_['lo']:.2f} with all reviews, {pp_['lo']:.2f} with the published deadlines and {p6_['lo']:.2f} without the 12-day cohort.", "lrrrrrrr")
ctn = {"2024-2026 (as in the calibration)": "2024--2026, pooled (as in the manuscript)", "2020-2022 (before the 2024 form)": "2020--2022, pooled", "2018-2022": "2018--2022, pooled", "2022 (adjacent, same scores)": "2022 (adjacent cohort)"}
lines = []
for k, lab in ctn.items():
    c_ = F4["contrasts"][k]
    lines.append(f"{lab} & {f(c_['gap'])} ({c_['se']:.3f}) & {c_['T_ref']:.1f} & {c_['eta']:.2f} \\\\")
lines.append("\\midrule")
for y_ in sorted(F4["single_cohorts"], key=int):
    c_ = F4["single_cohorts"][y_]
    lines.append(f"{y_} & {f(c_['gap'])} ({c_['se']:.3f}) & {c_['T']:.0f} & {c_['eta']:.2f} \\\\")
table("The compression elasticity $\\eta$ implied by the final-24-hour length gap of the 12-day cohort (2023) against other cohorts", "tab:fc-eta",
      "Reference cohorts & Length gap (s.e.) & Window (days) & $\\eta$", lines,
      f"Length gap: within-manuscript coefficient of the final 24 hours on log words (reference: four days or more before the deadline), pooled over cohorts with inverse-variance weights; for 2023 it is {f(F4['gap_2023'])} ({F4['se_2023']:.3f}). $\\eta=\\ln(\\text{{gap}}_{{2023}}/\\text{{gap}}_{{\\text{{ref}}}})/\\ln(T_{{\\text{{ref}}}}/12)$, with $T_{{\\text{{ref}}}}$ the mean window of the reference cohorts (for 2018--2022, of those with a documented window). The contrast with 2024--2026 crosses the change of review form and the arrival of language models in 2024, after which the length gap is smaller for reasons unrelated to the window (Table~\\ref{{tab:s14}}); the contrasts with the earlier cohorts do not. The value 1.3 used in the scenario table of the manuscript lies within the range of the first group and above that of the second.", "lrrr")
lines = [f"{r_['eta']:.2f} & {r_['0.0025']} & {r_['0.005']} & {r_['0.01']} \\\\" for r_ in F4["conference_optimum_by_eta"]]
table("Conference scenario: the window that maximises $V(T)$ as a function of the compression elasticity $\\eta$", "tab:fc-opt",
      "$\\eta$ & $T^*$, $\\delta$=0.0025 & $T^*$, $\\delta$=0.005 & $T^*$, $\\delta$=0.01", lines,
      f"Committed reviewers ($a\\equiv1$, $L=0$), delivery profile of 2024--2026, shortfalls of {F4['Delta_final48h']:.3f} for reports of the final 48 hours and {F4['Delta_late']:.3f} for late reports, $\\Delta(T)=\\Delta\\,(20.5/T)^{{\\eta}}$; $T^*$ searched on 3--90 days. The rows $\\eta=0$ and $\\eta=1.3$ are Panel~A of the scenario table of the manuscript.", "rrrr")
sc_ = F6["score"]; dv_ = F6["dev_panel"]; bn_ = sc_["2018-2026"]["bin_n"]; k1, k2 = "2018-2025 (1-10 scale)", "2026 (0-10 scale in steps of two)"
eq_ = R["A9"]["dev_panel_equivalence"]; n_ = lambda v: f"{v:,}".replace(",", "{,}")
other = ("\\paragraph{Checks reported in the text.} "
    f"\\emph{{The shortfall capped.}} The rule $\\Delta(T)=\\Delta\\,(20.5/T)^{{\\eta}}$ grows without bound as the window shortens: with $\\Delta=0.085$ and $\\eta=1.3$ it is {F7['Delta_T']['7']:.2f} at seven days and {F7['Delta_T']['3']:.2f}, more than the whole value of a report, at three; at the upper corner of the joint ranges ($\\Delta=0.30$, $\\eta=2$) it exceeds the value of a report for windows of 11 days or fewer. "
    f"Capping it at the value of a report changes none of the six conference optima of the scenario table and none of the {n_(F8['n'])} joint optima of either scenario, because the optimum never lies where the cap binds. The constant 20.5 is the mean window, in days, of the 2024--2026 cohorts on which the delivery profile and the reference shortfall are estimated. "
    f"\\emph{{The functional form in the joint analysis.}} The joint analysis of Tables~\\ref{{tab:s12}}--\\ref{{tab:s13}} varies the curvature of the exponential acceptance curve. When the functional form is drawn instead (each of the four curves of Table~\\ref{{tab:s7}} with probability one quarter, the other inputs as before), the journal optimum has a median of {F8['q']['50']:.0f} days, an interquartile range of {F8['q']['25']:.0f}--{F8['q']['75']:.0f} and a central 90\\% range of {F8['q']['5']:.0f}--{F8['q']['95']:.0f}; {100 * F8['share_at_90']:.1f}\\% of the draws are at the upper end of the search (90 days), nearly all of them under the straight line, and the form accounts on its own for {100 * F8['variance_share_form']:.0f}\\% of the variance. Medians by form (5th--95th percentiles): "
    + ", ".join(f"{nm_} {F8['by_form'][k_]['q']['50']:.0f} ({F8['by_form'][k_]['q']['5']:.0f}--{F8['by_form'][k_]['q']['95']:.0f})" for k_, nm_ in (("exponential", "exponential"), ("hyperbolic", "hyperbolic"), ("logarithmic", "logarithmic"), ("linear", "straight line"))) + " days. "
    f"\\emph{{Scores on two scales.}} The pooled final-24-hour coefficient on the overall score ({sc_['2018-2026']['coef']['last24h']['b']:.3f}) mixes two rating scales: it is {sc_[k1]['coef']['last24h']['b']:.3f} (s.e.\\ {sc_[k1]['coef']['last24h']['se']:.3f}) in 2018--2025, on the 1--10 scale, and {sc_[k2]['coef']['last24h']['b']:.3f} ({sc_[k2]['coef']['last24h']['se']:.3f}) in 2026, on the scale from 0 to 10 in steps of two. The corresponding coefficients on the absolute deviation from the co-reviewers' mean are {f(dv_[k1]['coef']['last24h']['b'])} ({dv_[k1]['coef']['last24h']['se']:.3f}) and {f(dv_[k2]['coef']['last24h']['b'])} ({dv_[k2]['coef']['last24h']['se']:.3f}). "
    f"\\emph{{Equivalence.}} The 90\\% interval of the pooled deviation coefficient is {f(eq_['ci90'][0])} to {eq_['ci90'][1]:.3f}; equivalence margins of 0.03, 0.05 and 0.10 points were examined, and the manuscript reports the smallest, about 2\\% of the early mean. "
    f"\\emph{{Reviews by timing.}} The within-manuscript regressions of the manuscript use {n_(sc_['2018-2026']['n'])} reviews: {n_(bn_['early'])} delivered four days or more before the deadline, {n_(bn_['mid'])} two to four days before, {n_(bn_['last48h'])} between 24 and 48 hours before, {n_(bn_['last24h'])} in the final 24 hours and {n_(bn_['late'])} late.\n")
final_sec = ("\\clearpage\n\\section*{Checks added in the final version}\n\\begin{sloppypar}\n\\noindent The tables of this section were added after a last internal review of the revised manuscript. They are produced by \\path{code/final_checks.py} from the archived outputs of the original pipeline and, for the published deadlines and the scores by scale, from the review-level file.\n\\end{sloppypar}\n\n"
             + other + "\nTables~\\ref{tab:fc-haz}--\\ref{tab:fc-opt} follow.\n\n" + "\n".join(P) + "\n")
P[:] = P_main
_tm = TAB.parent / "OnlineResource3_tmlr.tex"
_tmlr = ("\\clearpage\n" + _tm.read_text()) if _tm.exists() else ""
(TAB.parent / "OnlineResource3.tex").write_text(head + "\n".join(P) + _tmlr + final_sec + "\n\\end{document}\n")
print("tables written:", sorted(p.name for p in TAB.iterdir()), "| OnlineResource3.tex")
