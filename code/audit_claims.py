"""Audit of the revised manuscript and of the response letter.
usage: python code/audit_claims.py <repo> <output/revision> <manuscript folder> [<output/tmlr>]
The manuscript folder holds deadline_paper.tex (or deadline_paper_marked.tex), references.bib and tables/;
if ResponseToReviewers.tex is there too, the claims and the quotations of the letter are checked as well.
(1) every number introduced in the revision is recomputed from the JSON outputs and
    the sentence that states it is located in the manuscript;
(2) every passage that the letter quotes from the manuscript exists in it;
(3) citations and bibliography are in bijection; labels and references resolve."""
import json
import re
import sys
from pathlib import Path

REPO, REV, R1 = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
RES = json.loads((REPO / "output" / "results.json").read_text())
R = json.loads((REV / "revision_results.json").read_text())
_mk = R1 / "deadline_paper_marked.tex"
tex = (_mk if _mk.exists() else R1 / "deadline_paper.tex").read_text(encoding="utf8")
plain = re.sub(r"\\rev\{", "", tex)
_lt = R1 / "ResponseToReviewers.tex"
letter = _lt.read_text(encoding="utf8") if _lt.exists() else None   # the response letter is not in the public repository
fails = []
n = 0


def chk(desc, value, claim, text, tol=None, where=plain):
    """value: recomputed number; claim: number as printed; text: sentence fragment that must be in the manuscript."""
    global n
    if where is None:   # claim of the response letter, letter not present
        return
    n += 1
    ok_num = abs(value - claim) <= (tol if tol is not None else 0.5 * 10 ** (-len(str(claim).split(".")[1]) if "." in str(claim) else 0) + 1e-9)
    ok_txt = text in where
    if not (ok_num and ok_txt):
        fails.append((desc, value, claim, ok_num, ok_txt, text[:70]))


a1 = {r["stat"]: r for r in R["A1"]["table"]}
a7, a9, a12 = R["A7"], R["A9"], R["A12"]
# --- cohort-level inference
chk("median slope lo", a1["s_q50"]["lo"], -0.004, r"The 95\% interval of the median slope is $-0.004$ to 0.146")
chk("median slope hi", a1["s_q50"]["hi"], 0.146, r"$-0.004$ to 0.146 (exact permutation $p=0.056$)")
chk("perm p", a1["s_q50"]["p_perm"], 0.056, r"exact permutation $p=0.056$")
chk("q25 lo", a1["s_q25"]["lo"], -0.022, r"that of the lower quartile $-0.022$ to 0.048")
chk("q25 hi", a1["s_q25"]["hi"], 0.048, r"$-0.022$ to 0.048")
chk("q90 lo", a1["s_q90"]["lo"], 0.072, r"that of the ninetieth percentile 0.072 to 0.717")
chk("q90 hi", a1["s_q90"]["hi"], 0.717, r"0.072 to 0.717")
chk("passthrough lower bound", 1 - a1["s_q50"]["hi"], 0.85, r"a pass-through of the window to the delivery date below 0.85")
chk("atom slope pp", 100 * a1["share_atom"]["slope"], -0.6, r"($-0.6$ percentage points per day of window; interval $-1.6$ to 0.4)")
chk("atom lo pp", 100 * a1["share_atom"]["lo"], -1.6, r"interval $-1.6$ to 0.4")
chk("atom hi pp", 100 * a1["share_atom"]["hi"], 0.4, r"interval $-1.6$ to 0.4")
chk("late slope pp", 100 * a1["share_late"]["slope"], -1.2, r"of the late share ($-1.2$; $-3.3$ to 1.0)")
chk("late lo pp", 100 * a1["share_late"]["lo"], -3.3, r"($-1.2$; $-3.3$ to 1.0)")
chk("late hi pp", 100 * a1["share_late"]["hi"], 1.0, r"($-1.2$; $-3.3$ to 1.0)")
chk("leverage 2023", R["A1"]["leverage"]["2023"], 0.77, r"has a leverage of 0.77")
chk("loo min", a1["s_q50"]["loo_min"], 0.06, r"leaves the median slope between 0.06 and 0.10")
chk("loo max", a1["s_q50"]["loo_max"], 0.10, r"between 0.06 and 0.10")
chk("slope since assignment hi (interval 0.85 to 1.00)", 1 - a1["s_q50"]["lo"], 1.00, r"by 0.93 days per day at the median of on-time reviews with an interval of 0.85 to 1.00")
chk("rct pass-through", R["A4"]["rho_rct"], 0.71, r"cut the median turnaround by 12 days, a pass-through of 0.71")
e3 = {r["year"]: r for r in RES["E3"]["table"]}
chk("scale factor ~10", e3[2026]["n"] / e3[2020]["n"], 10, r"by a factor of about ten in the number of reviews", tol=0.5)
# --- equivalence, before/after, kappa
eq = a9["dev_panel_equivalence"]
chk("TOST p", eq["tost"]["0.03"]["p"], 0.010, r"(two one-sided tests, $p=0.010$;")
chk("margin pct", eq["tost"]["0.03"]["margin_pct_of_mean"], 2, r"about 2\% of the early mean of 1.37", tol=0.5)
chk("early mean", eq["mean_early"], 1.37, r"early mean of 1.37")
chk("dev coef", eq["coef"], -0.010, r"is 0.010 lower than an early review's, with a standard error of 0.008")
chk("pre length", a9["log_words"]["pre"]["est"], -0.119, r"is $-0.119$ (standard error 0.010) in 2018--2023")
chk("pre length se", a9["log_words"]["pre"]["se"], 0.010, r"$-0.119$ (standard error 0.010)")
chk("post length", a9["log_words"]["post"]["est"], -0.078, r"and $-0.078$ (0.004) in 2024--2026")
chk("post length se", a9["log_words"]["post"]["se"], 0.004, r"$-0.078$ (0.004)")
chk("pre conf", a9["confidence"]["pre"]["est"], -0.099, r"($-0.099$ and $-0.129$; Online Resource~3, Table~S16)")
chk("post conf", a9["confidence"]["post"]["est"], -0.129, r"($-0.099$ and $-0.129$;")
chk("post excl 2025", a9["log_words"]["post_excl_2025"]["est"], -0.078, r"leaves the estimate unchanged ($-0.078$)")
pl = R["A10"]["pooled"]
chk("kappa EE-LL", pl["early-early minus late-late | 2021-2025"]["est"], 0.06, r"pooled over these five it is 0.06 (0.02 to 0.10)")
chk("kappa EE-AA", pl["early-early minus atom-atom | 2021-2025"]["est"], 0.03, r"final-48-hour pairs is 0.03")
KB = json.loads((REV / "kappa_bootstrap.json").read_text())
assert KB["max_abs_diff_with_archived_kappa"] < 1e-12
kb1, kb2 = KB["pooled"]["early-early minus late-late | 2021-2025"], KB["pooled"]["early-early minus atom-atom | 2021-2025"]
chk("kappa EE-LL bootstrap weights", kb1["est"], 0.06, r"pooled over these five it is 0.06 (0.02 to 0.10)")
chk("kappa EE-LL lo", kb1["lo"], 0.02, r"pooled over these five it is 0.06 (0.02 to 0.10)")
chk("kappa EE-LL hi", kb1["hi"], 0.10, r"pooled over these five it is 0.06 (0.02 to 0.10)")
kb8 = KB["pooled"]["early-early minus late-late | 2018-2025"]
chk("kappa all cohorts", kb8["est"], 0.05, r"late--late pairs agree less than early--early pairs by 0.05, with a 95\% interval of 0.02 to 0.09 from resampling manuscripts")
chk("kappa all cohorts lo", kb8["lo"], 0.02, r"a 95\% interval of 0.02 to 0.09 from resampling manuscripts")
chk("kappa all cohorts hi", kb8["hi"], 0.09, r"a 95\% interval of 0.02 to 0.09 from resampling manuscripts")
chk("kappa all cohorts (intro)", kb8["est"], 0.05, r"(Cohen's $\kappa$ lower by about 0.05)")
chk("kappa EE-AA bootstrap weights", kb2["est"], 0.03, r"final-48-hour pairs is 0.03 (0.003 to 0.057;")
chk("kappa EE-AA lo", kb2["lo"], 0.003, r"final-48-hour pairs is 0.03 (0.003 to 0.057;")
chk("kappa EE-AA hi", kb2["hi"], 0.057, r"final-48-hour pairs is 0.03 (0.003 to 0.057;")
chk("letter kappa all", kb8["est"], 0.05, r"late--late pairs agree less than early--early pairs by 0.05, with a 95\% interval of 0.02 to 0.09 from a bootstrap over manuscripts", where=letter)
chk("letter kappa lo", kb8["lo"], 0.02, r"a 95\% interval of 0.02 to 0.09 from a bootstrap over manuscripts", where=letter)
chk("letter kappa hi", kb8["hi"], 0.09, r"a 95\% interval of 0.02 to 0.09 from a bootstrap over manuscripts", where=letter)
by = pl["early-early minus late-late | 2021-2025"]["by_year"]
chk("kappa min diff (2022)", by["2022"] if "2022" in by else by[2022], 0.005, r"though by only 0.005 in 2022 and by as much as 0.13 in 2023")
chk("kappa min is 2022", min(by.values()), 0.005, r"though by only 0.005 in 2022")
chk("kappa max diff (2023)", by["2023"] if "2023" in by else by[2023], 0.13, r"by as much as 0.13 in 2023")
chk("kappa max is 2023", max(by.values()), 0.13, r"by as much as 0.13 in 2023")
chk("kappa positive in 5", pl["early-early minus late-late | 2021-2025"]["n_positive"], 5, r"In the five cohorts with at least 225 early--early pairs (2021--2025) the difference is positive in each")
# --- calibration: scenario table
ts = RES["E7"]["tstar"]
jr = [r for r in ts if not r["conference"]]; cf = [r for r in ts if r["conference"]]
chk("journal min", min(r["T_star"] for r in jr), 23, r"The optimum is 23--44 days across the panel")
chk("journal max", max(r["T_star"] for r in jr), 44, r"23--44 days across the panel")
cen = sorted(r["T_star"] for r in jr if r["L"] == 0.3 and r["delta"] == 0.005)
chk("journal central lo", cen[0], 32, r"(32--35 at the central values $L=0.3$ and $\delta=0.005$)")
chk("journal central hi", cen[1], 35, r"(32--35 at the central values")
c13 = sorted(r["T_star"] for r in cf if r["eta"] == 1.3)
chk("conf eta1.3 lo", c13[0], 14, r"the optimum is interior at 14--25 days (19 at the central delay cost)")
chk("conf eta1.3 hi", c13[-1], 25, r"interior at 14--25 days")
chk("conf central", c13[1], 19, r"(19 at the central delay cost)")
c0 = sorted(r["T_star"] for r in cf if r["eta"] == 0.0)
chk("conf eta0 range", c0[0] * 10 + c0[-1], 37, r"(3--7 days with our delay costs)")
chk("loss28 conf eta0 lo", a12["loss28_confTrue_eta0.0"][0], 4, r"With $\eta=0$ a 28-day window costs 4--20\% of a report", tol=0.5)
chk("loss28 conf eta0 hi", a12["loss28_confTrue_eta0.0"][1], 20, r"a 28-day window costs 4--20\%", tol=0.5)
chk("loss14/21 conf eta1.3 max", max(a12["loss14_confTrue_eta1.3"][1], a12["loss21_confTrue_eta1.3"][1]), 2.5, r"the loss from 14 or 21 days is at most 2.5\%")
chk("loss42 conf eta1.3 lo", a12["loss42_confTrue_eta1.3"][0], 2, r"the loss from 42 days is 2--19\%", tol=0.5)
chk("loss42 conf eta1.3 hi", a12["loss42_confTrue_eta1.3"][1], 19, r"2--19\%", tol=0.5)
chk("loss7 conf eta1.3 lo", a12["loss7_confTrue_eta1.3"][0], 9, r"a 7-day window costs 9--17\% of a report and a 42-day window 2--19\%", tol=0.5)
chk("loss7 conf eta1.3 hi", a12["loss7_confTrue_eta1.3"][1], 17, r"9--17\%", tol=0.5)
l14 = [a12["loss14_confFalse_eta0.0"][0], a12["loss14_confFalse_eta1.3"][1]]
chk("loss14 journal lo", l14[0], 4, r"a 14-day window costs 4--24\% of a report", tol=0.5)
chk("loss14 journal hi", l14[1], 24, r"4--24\%", tol=0.5)
chk("loss28 journal max", max(a12["loss28_confFalse_eta0.0"][1], a12["loss28_confFalse_eta1.3"][1]), 3.6, r"a 28-day window at most 3.6\%", tol=0.05)
chk("loss42 journal max", max(a12["loss42_confFalse_eta0.0"][1], a12["loss42_confFalse_eta1.3"][1]), 7, r"and a 42-day window at most 7\%", tol=0.5)
r32 = [r for r in jr if r["L"] == 0.3 and r["delta"] == 0.005 and r["eta"] == 0.0][0]
abar, lam = R["A2"]["abar"], R["A2"]["lam_a"]
import math
a14 = abar * (1 - math.exp(-lam * 14))
chk("declined share at 14 d", (r32["a_Tstar"] - a14) / r32["a_Tstar"], 0.2, r"about a fifth of the invitations that would have been accepted at the optimum are declined", tol=0.03)
chk("marginal gain point", a12["marginal_quality_per_day_pct"]["point"], 0.15, r"reduces the shortfall by 0.15\% of a report, less than any delay cost we consider")
chk("marginal gain upper", a12["marginal_quality_per_day_pct"]["upper"], 0.42, r"the gain would be 0.42\% of a report per day, more than our lowest delay cost")
chk("closed form", a12["closed_form_lam0"], 32.0, r"32.0}$ days", tol=0.05)
# --- alternatives
a3 = {(r["form"], r["delta"], r["eta"]): r["T_star"] for r in R["A3"]["table"]}
forms = ["exponential (submitted)", "hyperbolic", "logarithmic", "linear (chord)"]
for d, claim, txt in ((0.005, [32, 30, 28, 13], r"the four curves give optima of 32, 30, 28 and 13 days"), (0.01, [25, 21, 16, 5], r"at $\delta=0.01$, 25, 21, 16 and 5 days"), (0.0025, [40, 44, 51, 90], r"at $\delta=0.0025$, 40, 44, 51 and 90 days or more")):
    for f_, c in zip(forms, claim):
        chk(f"acceptance form {f_} delta {d}", a3[(f_, d, 0.0)], c, txt)
lamt = {(round(r["lam_a"], 3), r["delta"]): r["T_star"] for r in R["A3_lambda"]}
chk("lam 0.05 lo", lamt[(0.05, 0.01)], 35, r"$\lambda_a=0.05$ gives 35--63 days"); chk("lam 0.05 hi", lamt[(0.05, 0.0025)], 63, r"35--63 days")
chk("lam 0.20 lo", lamt[(0.2, 0.01)], 16, r"$\lambda_a=0.20$ gives 16--24"); chk("lam 0.20 hi", lamt[(0.2, 0.0025)], 24, r"gives 16--24")
pc = R["A4"]["closed"]
r71 = [r["T_closed"] for r in pc if abs(r["rho"] - R["A4"]["rho_rct"]) < 1e-9]; r1_ = [r["T_closed"] for r in pc if r["rho"] == 1.0]
chk("rct closed lo", min(r71), 25, r"gives 25--45 days, against 22--41 under full anchoring", tol=0.5); chk("rct closed hi", max(r71), 45, r"25--45 days", tol=0.5)
chk("full anchoring lo", min(r1_), 22, r"against 22--41", tol=0.5); chk("full anchoring hi", max(r1_), 41, r"22--41", tol=0.5)
a5 = {r["tail"]: r for r in R["A5"]["table"]}
cols = [c for c in a5["power"] if c.startswith("T_")]
d_exp = max(abs(a5["exponential"][c] - a5["power"][c]) for c in cols); d_tr = max(abs(a5["truncated"][c] - a5["power"][c]) for c in cols); d_fl = max(abs(a5["floor"][c] - a5["power"][c]) for c in cols)
chk("tails exp/trunc max shift", max(d_exp, d_tr), 1, r"moves the optimum by at most one day")
chk("tail floor max shift", d_fl, 3, r"moves it by up to three days")
a6 = {(r["Delta"], r["eta"], r["delta"]): r for r in R["A6"]}
chk("Delta0 conf", max(a6[(0.0, e, d)]["T_conference"] for e in (0.0, 1.3) for d in (0.0025, 0.005, 0.01)), 3, r"the conference optimum is the shortest window considered (3 days)")
j0 = [a6[(0.0, e, d)]["T_journal"] for e in (0.0, 1.3) for d in (0.0025, 0.005, 0.01)]
chk("Delta0 journal lo", min(j0), 25, r"and the journal optimum at $L=0.3$ is 25--40 days, as with $\Delta=0.085$"); chk("Delta0 journal hi", max(j0), 40, r"25--40 days")
j85 = [a6[(0.085, 0.0, d)]["T_journal"] for d in (0.0025, 0.005, 0.01)]
chk("Delta .085 journal lo", min(j85), 25, r"is 25--40 days, as with $\Delta=0.085$"); chk("Delta .085 journal hi", max(j85), 40, r"is 25--40 days, as with $\Delta=0.085$")
chk("journal rise with Delta, eta 1.3", max(a6[(0.3, 1.3, d)]["T_journal"] - a6[(0.0, 1.3, d)]["T_journal"] for d in (0.0025, 0.005, 0.01)), 9, r"and rises by up to nine days when $\eta=1.3$")
for d in (0.0025, 0.005, 0.01):
    chk(f"journal fall with Delta, eta 0, delta {d}", a6[(0.0, 0.0, d)]["T_journal"] - a6[(0.3, 0.0, d)]["T_journal"], 1, r"the journal optimum falls by a day when $\eta=0$")
c3 = [a6[(0.3, 1.3, d)]["T_conference"] for d in (0.0025, 0.005, 0.01)]
chk("conf Delta .3 lo", min(c3), 23, r"rises from 3 days to 23--40 when $\eta=1.3$"); chk("conf Delta .3 hi", max(c3), 40, r"23--40 when")
# --- joint
cq, jq = a7["conference"]["q"], a7["journal"]["q"]
chk("joint n", a7["n"], 20000, r"We drew 20,000 combinations of the inputs")
chk("conf median", cq["50"], 18, r"For the conference the optimum has a median of 18 days, an interquartile range of 11--25 and a central 90\% range of 5--33")
chk("conf q25", cq["25"], 11, r"interquartile range of 11--25"); chk("conf q75", cq["75"], 25, r"11--25 and a central")
chk("conf q5", cq["5"], 5, r"central 90\% range of 5--33"); chk("conf q95", cq["95"], 33, r"range of 5--33")
chk("conf eta<0.5 median", a7["conference"]["eta<0.5"]["50"], 9, r"its median is 9 days when $\eta<0.5$ and 23 when $\eta\ge1$")
chk("conf eta>=1 median", a7["conference"]["eta>=1"]["50"], 23, r"and 23 when $\eta\ge1$")
chk("S conf Delta", 100 * a7["S_conference"]["Delta"], 41, r"$\Delta$ and $\eta$ account for 41\% and 39\% of its variance, the delay cost for 13\%", tol=0.5)
chk("S conf eta", 100 * a7["S_conference"]["eta"], 39, r"41\% and 39\%", tol=0.5); chk("S conf delta", 100 * a7["S_conference"]["delta"], 13, r"the delay cost for 13\%", tol=0.5)
chk("jour median", jq["50"], 37, r"For the journal the median is 37 days, the interquartile range 30--46 and the central 90\% range 23--60")
chk("jour q25", jq["25"], 30, r"30--46"); chk("jour q75", jq["75"], 46, r"30--46"); chk("jour q5", jq["5"], 23, r"23--60"); chk("jour q95", jq["95"], 60, r"23--60")
chk("S jour lam", 100 * a7["S_journal"]["lam_a"], 64, r"accounts for 64\% of the variance and the delay cost for 23\%", tol=0.5)
chk("S jour delta", 100 * a7["S_journal"]["delta"], 23, r"the delay cost for 23\%", tol=0.5)
chk("S jour others < 3", 100 * max(v for k, v in a7["S_journal"].items() if k not in ("lam_a", "delta", "interactions")), 2.5, r"and no other input for more than 3\%", tol=0.5)
chk("jour loss14>loss42", 100 * a7["journal"]["loss14_gt_loss42"], 95, r"In 95\% of the journal draws with exponential curves a 14-day window costs more than a 42-day one", tol=0.5)
chk("jour median loss14", a7["journal"]["loss14_q"]["50"], 17, r"(median losses of 17\% and 0.8\% of a report)", tol=0.5)
chk("jour median loss42", a7["journal"]["loss42_q"]["50"], 0.8, r"17\% and 0.8\%")
chk("conf loss14<loss42", 100 * (1 - a7["conference"]["loss14_gt_loss42"]), 71, r"in 71\% of the conference draws it costs less", tol=0.5)
chk("unimodal conf", a7["conference"]["unimodal"], 1.0, r"has a single peak on 3--90 days in every draw of both scenarios", tol=0); chk("unimodal jour", a7["journal"]["unimodal"], 1.0, r"in every draw of both scenarios", tol=0)
# --- data facts restated in new text
e1 = {r["year"]: r for r in RES["E1"]}
late = {y: 100 * e1[y]["share_late"] for y in e1}
chk("late 2018-24 min", min(late[y] for y in range(2018, 2025)), 21, r"the late share fell from 21--37\% in 2018--2024 to 16--18\% in 2025--2026", tol=0.5)
chk("late 2018-24 max", max(late[y] for y in range(2018, 2025)), 37, r"21--37\%", tol=0.5)
chk("late 2025", late[2025], 16, r"16--18\% in 2025--2026", tol=0.5); chk("late 2026", late[2026], 18, r"16--18\%", tol=0.5)
chk("late 2024", late[2024], 27, r"was no more punctual than its predecessors (27\% late)", tol=0.5)
med = {y: e1[y]["median_days_before"] for y in e1}
chk("median early min", min(med[y] for y in range(2018, 2025)), 0.5, r"the median review moved from 0.5--0.9 to 1.2--1.4 days before the deadline", tol=0.05)
chk("median early max", max(med[y] for y in range(2018, 2025)), 0.9, r"0.5--0.9", tol=0.05)
chk("median 2026", med[2026], 1.2, r"1.2--1.4", tol=0.05); chk("median 2025", med[2025], 1.4, r"1.2--1.4", tol=0.05)
mu = {r["year"]: r["mu"] for r in RES["E2"]} if "mu" in RES["E2"][0] else None
if mu:
    chk("mu last two min", min(mu[2025], mu[2026]), 0.38, r"0.38--0.42 in the last two cohorts against 0.20--0.39 in 2020--2024", tol=0.005)
    chk("mu last two max", max(mu[2025], mu[2026]), 0.42, r"0.38--0.42", tol=0.005)
    chk("mu 2020-24 min", min(mu[y] for y in range(2020, 2025)), 0.20, r"against 0.20--0.39 in 2020--2024", tol=0.005)
    chk("mu 2020-24 max", max(mu[y] for y in range(2020, 2025)), 0.39, r"0.20--0.39", tol=0.005)
else:
    fails.append(("mu not found in E2", list(RES["E2"][0].keys())))
chk("reviews 2026/2023", e1[2026]["n_reviews"] / e1[2023]["n_reviews"], 4, r"with four times as many reviews in 2026 as in 2023", tol=0.2)
a11 = {r["year"]: r for r in R["A11"]}
chk(">=3 reviews 2021-26", 100 * (1 - max(a11[y]["share_lt3"] for y in range(2021, 2027))), 99.8, r"(99.8\% or more in each cohort from 2021 to 2026;", tol=0.05)
# --- checks of the final version (final_checks.json)
FC = json.loads((REV / "final_checks.json").read_text())
F1, F2, F3, F4, F5, F6, F7, F8 = (FC[k] for k in ("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8"))
sl = F1["pooled_log_hazard_slope"]
chk("F1 hazard slope", sl["slope"], 0.055, r"the slope of the log hazard on the window is 0.055 (95\% interval $-0.027$ to 0.137; Table~S26)")
chk("F1 lo", sl["lo"], -0.027, r"(95\% interval $-0.027$ to 0.137; Table~S26)"); chk("F1 hi", sl["hi"], 0.137, r"(95\% interval $-0.027$ to 0.137; Table~S26)")
chk("gamma slope", a1["gamma"]["slope"], -0.006, r"the slope of the anchoring rate $\hat\gamma$ of the mixture is $-0.006$ ($-0.016$ to 0.003)")
chk("gamma lo", a1["gamma"]["lo"], -0.016, r"is $-0.006$ ($-0.016$ to 0.003)"); chk("gamma hi", a1["gamma"]["hi"], 0.003, r"is $-0.006$ ($-0.016$ to 0.003)")
for k, c, d in (("q25", 0.005, None), ("q50", 0.02, None), ("q75", 0.09, None), ("q90", 0.29, None)):
    chk(f"F2 implied {k}", F2["power"][k]["implied_slope"], c, r"gives slopes of 0.005, 0.02, 0.09 and 0.29 for the lower quartile, the median, the upper quartile and the ninetieth percentile")
for k, c in (("s_q25", 0.013), ("s_q50", 0.07), ("s_q75", 0.20), ("s_q90", 0.40)):
    chk(f"observed slope {k}", a1[k]["slope"], c, r"against the observed 0.013, 0.07, 0.20 and 0.40 (Online Resource~3, Table~S27)")
chk("slope q50", a1["s_q50"]["slope"], 0.071, r"the median number of days before the deadline moves 0.071 days per day of window, the lower quartile 0.013, and the median number of days since assignment 0.929")
chk("pass-through on-time", 1 - a1["s_q50"]["slope"], 0.929, r"and the median number of days since assignment 0.929; with late reviews included the last figure is 0.942")
pa = F3["all_reviews"]["q50"]
chk("F3 all reviews", pa["passthrough"], 0.942, r"with late reviews included the last figure is 0.942")
chk("F3 all reviews (intro)", pa["passthrough"], 0.94, r"0.93 days later after assignment (0.94 with late reviews included)")
chk("F3 all reviews letter", pa["passthrough"], 0.94, r"with late reviews included the pass-through of the median is 0.94 (0.87 to 1.01)", where=letter)
chk("F3 all lo letter", pa["lo"], 0.87, r"is 0.94 (0.87 to 1.01)", where=letter); chk("F3 all hi letter", pa["hi"], 1.01, r"is 0.94 (0.87 to 1.01)", where=letter)
p6 = F3["on_time_without_2023"]
chk("F3 without 2023 lo", p6["lo"], 0.73, r"that bound rests on the 12-day cohort (without it the interval is 0.73 to 1.07)")
chk("F3 without 2023 hi", p6["hi"], 1.07, r"(without it the interval is 0.73 to 1.07)")
pp = F5["published"]["passthrough_median"]
chk("F5 published", pp["est"], 0.89, r"is 0.89 (95\% interval 0.77 to 1.01) instead of 0.93 (0.85 to 1.00) (Online Resource~3, Table~S27)")
chk("F5 lo", pp["lo"], 0.77, r"is 0.89 (95\% interval 0.77 to 1.01)"); chk("F5 hi", pp["hi"], 1.01, r"is 0.89 (95\% interval 0.77 to 1.01)")
chk("F5 text", pp["lo"], 0.77, r"and on the deadlines we use (0.77 to 1.01 with the published ones)")
chk("F5 letter", pp["est"], 0.89, r"the pass-through of the median is 0.89 (0.77 to 1.01) instead of 0.93 (0.85 to 1.00)", where=letter)
chk("nine days later", 9 * (1 - a1["s_q50"]["slope"]), 8.4, r"go with a median delivery 8.4 days later after assignment and 0.6 days earlier before the deadline", tol=0.05)
chk("nine days earlier", 9 * a1["s_q50"]["slope"], 0.6, r"8.4 days later after assignment and 0.6 days earlier before the deadline", tol=0.05)
c_ = F4["contrasts"]; sg = F4["single_cohorts"]
chk("eta later", c_["2024-2026 (as in the calibration)"]["eta"], 1.2, r"an elasticity $\eta$ of about 1.2 with respect to the ratio of the windows (1.0 to 1.6 for the three cohorts taken singly)", tol=0.05)
chk("eta later min", min(sg[y]["eta"] for y in ("2024", "2025", "2026")), 1.0, r"(1.0 to 1.6 for the three cohorts taken singly)", tol=0.05)
chk("eta later max", max(sg[y]["eta"] for y in ("2024", "2025", "2026")), 1.6, r"(1.0 to 1.6 for the three cohorts taken singly)", tol=0.05)
chk("eta earlier", c_["2020-2022 (before the 2024 form)"]["eta"], 0.7, r"the elasticity is about 0.7 (0.5 to 0.9 singly)", tol=0.05)
chk("eta earlier min", min(sg[y]["eta"] for y in ("2020", "2021", "2022")), 0.5, r"about 0.7 (0.5 to 0.9 singly)", tol=0.05)
chk("eta earlier max", max(sg[y]["eta"] for y in ("2020", "2021", "2022")), 0.9, r"about 0.7 (0.5 to 0.9 singly)", tol=0.05)
chk("gap 2023", F4["gap_2023"], -0.146, r"the 2023 gap ($-0.146$) implies")
chk("gap 2020-22", c_["2020-2022 (before the 2024 form)"]["gap"], -0.107, r"2020--2022 ($-0.107$ pooled)")
opt = {r_["eta"]: r_ for r_ in F4["conference_optimum_by_eta"]}
chk("opt eta .7 lo", opt[0.7]["0.01"], 9, r"With the elasticity suggested by the earlier cohorts ($\eta=0.7$) the optimum is 9--19 days (14 at the central delay cost)")
chk("opt eta .7 hi", opt[0.7]["0.0025"], 19, r"the optimum is 9--19 days (14 at the central delay cost)")
chk("opt eta .7 central", opt[0.7]["0.005"], 14, r"(14 at the central delay cost)")
chk("opt eta 0 central", opt[0.0]["0.005"], 5, r"from 5 days at $\eta=0$ to 22 at $\eta=2$ (Online Resource~3, Table~S29)")
chk("opt eta 2 central", opt[2.0]["0.005"], 22, r"from 5 days at $\eta=0$ to 22 at $\eta=2$")
chk("opt intro .7", opt[0.7]["0.01"] * 100 + opt[0.7]["0.0025"], 919, r"($\eta\approx0.7$), the optimum is 9--19 days; under the comparison with the later cohorts ($\eta\approx1.3$) it is 14--25 days")
chk("opt two conclusions", opt[1.3]["0.005"], 19, r"5 days if last-minute losses do not grow as windows shorten, 14 if they grow as the earlier cohorts suggest, and 19 if they double between 21 and 12 days")
chk("opt letter", opt[2.0]["0.005"], 22, r"5 days for $\eta=0$, 14 for $\eta=0.7$, 19 for $\eta=1.3$ and 22 for $\eta=2$", where=letter)
sc = F6["score"]; dv = F6["dev_panel"]; k1, k2 = "2018-2025 (1-10 scale)", "2026 (0-10 scale in steps of two)"
chk("score 2018-25", sc[k1]["coef"]["last24h"]["b"], 0.13, r"(0.13 on the 1--10 scale of 2018--2025 and 0.22 on the 2026 scale")
chk("score 2026", sc[k2]["coef"]["last24h"]["b"], 0.22, r"0.13 on the 1--10 scale of 2018--2025 and 0.22 on the 2026 scale")
chk("score intro", sc[k1]["coef"]["last24h"]["b"], 0.13, r"give higher scores (by 0.13 points on the 1--10 scale used until 2025)")
T6 = (R1 / "tables" / "tab4_quality.tex").read_text()
bn = sc["2018-2026"]["bin_n"]
for b_, c in (("early", 43269), ("mid", 29595), ("last48h", 32018), ("last24h", 58536), ("late", 45959)):
    chk(f"bin n {b_}", bn[b_], c, r"(43,269 early, 29,595 two to four days before, 32,018 between 24 and 48 hours before, 58,536 in the final 24 hours and 45,959 late)", where=T6)
chk("T6 score scale 1", sc[k1]["coef"]["last24h"]["b"], 0.130, r"are 0.130 (0.015) and $-0.017$ (0.010)", where=T6); chk("T6 score se", sc[k1]["coef"]["last24h"]["se"], 0.015, r"are 0.130 (0.015)", where=T6)
chk("T6 dev scale 1", dv[k1]["coef"]["last24h"]["b"], -0.017, r"0.130 (0.015) and $-0.017$ (0.010)", where=T6); chk("T6 dev se", dv[k1]["coef"]["last24h"]["se"], 0.010, r"$-0.017$ (0.010)", where=T6)
chk("T6 score scale 2", sc[k2]["coef"]["last24h"]["b"], 0.224, r"0.224 (0.023) and $-0.001$ (0.015)", where=T6); chk("T6 score se 2", sc[k2]["coef"]["last24h"]["se"], 0.023, r"0.224 (0.023)", where=T6)
chk("T6 dev scale 2", dv[k2]["coef"]["last24h"]["b"], -0.001, r"$-0.001$ (0.015)", where=T6); chk("T6 dev se 2", dv[k2]["coef"]["last24h"]["se"], 0.015, r"$-0.001$ (0.015)", where=T6)
chk("dev 90 lo", eq["ci90"][0], -0.024, r"the 90\% interval of the difference is $-0.024$ to 0.003"); chk("dev 90 hi", eq["ci90"][1], 0.003, r"$-0.024$ to 0.003")
chk("dev pre-2024", a9["dev_panel"]["pre"]["est"], -0.043, r"Before 2024 last-minute reviews deviated less than early ones ($-0.043$, standard error 0.016; Online Resource~3, Table~S16)")
chk("dev pre-2024 se", a9["dev_panel"]["pre"]["se"], 0.016, r"($-0.043$, standard error 0.016;")
chk("F8 median", F8["q"]["50"], 39, r"the median is 39 days, the interquartile range 32--51 and the central 90\% range runs from 21 days to the 90 that bound the search")
chk("F8 q25", F8["q"]["25"], 32, r"the interquartile range 32--51"); chk("F8 q75", F8["q"]["75"], 51, r"the interquartile range 32--51")
chk("F8 q5", F8["q"]["5"], 21, r"runs from 21 days to the 90 that bound the search"); chk("F8 q95", F8["q"]["95"], 90, r"runs from 21 days to the 90 that bound the search")
chk("F8 letter", F8["q"]["5"] * 100 + F8["q"]["95"], 2190, r"the central 90\% range widens to 21--90 days", where=letter)
chk("F7 cap unchanged conf", F7["conference"]["share_changed"], 0.0, r"Capping it at that value changes none of the optima we report (Online Resource~3)", tol=0)
chk("F7 cap unchanged journal", F7["journal"]["share_changed"], 0.0, r"Capping it at that value changes none of the optima we report", tol=0)
chk("F7 cap table", sum(abs(r_["capped"] - r_["uncapped"]) for r_ in F7["table10_conference"]), 0, r"Capping it at that value changes none of the optima we report", tol=0)
chk("F7 Delta(3) > 1", F7["Delta_T"]["3"], 1.04, r"at three days it exceeds the whole value of a report", tol=0.01)
chk("ref window 20.5", sum(e3[y]["T"] for y in (2024, 2025, 2026)) / 3, 20.5, r"where 20.5 days is the mean window of the cohorts on which $\Delta$ is measured", tol=0.2)
chk("Delta final 48 h", RES["E7"]["Delta_atom"], 0.085, r"that of the final 24 hours, 0.085, for reports delivered in the final 48 hours, and that of late reports, 0.095, for late reports")
chk("Delta late", RES["E7"]["Delta_late"], 0.095, r"and that of late reports, 0.095, for late reports")
val = {int(v["year"]): v for v in RES["E7"]["validation"]}; oos = [val[y] for y in (2020, 2021, 2022, 2023)]
gap = [100 * (v["pred_share_pre"] - v["obs_share_pre"]) for v in oos]
chk("oos early gap min", min(gap), 9, r"is over-predicted by 9 to 16 points (36--39\% against 20--30\% observed)", tol=0.6); chk("oos early gap max", max(gap), 16, r"over-predicted by 9 to 16 points", tol=0.6)
chk("oos early pred min", 100 * min(v["pred_share_pre"] for v in oos), 36, r"(36--39\% against 20--30\% observed)", tol=0.6); chk("oos early pred max", 100 * max(v["pred_share_pre"] for v in oos), 39, r"(36--39\% against 20--30\% observed)", tol=0.6)
chk("oos early obs min", 100 * min(v["obs_share_pre"] for v in oos), 20, r"(36--39\% against 20--30\% observed)", tol=0.6); chk("oos early obs max", 100 * max(v["obs_share_pre"] for v in oos), 30, r"(36--39\% against 20--30\% observed)", tol=0.6)
chk("oos late pred min", 100 * min(v["pred_share_late"] for v in oos), 19, r"the late share is under-predicted (19--20\% against 21--34\%)", tol=0.6); chk("oos late pred max", 100 * max(v["pred_share_late"] for v in oos), 20, r"(19--20\% against 21--34\%)", tol=0.6)
chk("oos late obs min", 100 * min(v["obs_share_late"] for v in oos), 21, r"(19--20\% against 21--34\%)", tol=0.6); chk("oos late obs max", 100 * max(v["obs_share_late"] for v in oos), 34, r"(19--20\% against 21--34\%)", tol=0.6)
chk("oos 2023 late", 100 * val[2023]["obs_share_late"], 34, r"most of all for the 12-day cohort (20\% against 34\%; Online Resource~1, Table~S4)", tol=0.6); chk("oos 2023 late pred", 100 * val[2023]["pred_share_late"], 20, r"(20\% against 34\%;", tol=0.6)
chk("oos atom pred min", 100 * min(v["pred_share_atom"] for v in oos), 42, r"(predicted 42--44\%, observed 43--50\%)", tol=0.6); chk("oos atom pred max", 100 * max(v["pred_share_atom"] for v in oos), 44, r"(predicted 42--44\%, observed 43--50\%)", tol=0.6)
chk("oos atom obs min", 100 * min(v["obs_share_atom"] for v in oos), 43, r"(predicted 42--44\%, observed 43--50\%)", tol=0.6); chk("oos atom obs max", 100 * max(v["obs_share_atom"] for v in oos), 50, r"(predicted 42--44\%, observed 43--50\%)", tol=0.6)
chk("growth submissions", e1[2026]["n_papers"] / e1[2018]["n_papers"], 20, r"twenty-fold growth in submissions (twenty-five-fold in reviews) and five changes of review form", tol=0.5)
chk("growth reviews", e1[2026]["n_reviews"] / e1[2018]["n_reviews"], 25.6, r"(twenty-five-fold in reviews)", tol=0.1)
chk("a(14)", abar * (1 - math.exp(-lam * 14)), 0.51, r"which falls steeply below four weeks ($a(14)=0.51$)")
chk("hazard at T*=32", 100 * 1.10 * 33 ** -1.85, 0.17, r"the hazard is about 0.17\% per day and $a'/a$ about 0.44\%")
chk("a'/a at 32", 100 * lam * math.exp(-lam * 32) / (1 - math.exp(-lam * 32)), 0.44, r"$a'/a$ about 0.44\%")
print(f"(1) numeric claims checked: {n}; failures: {len(fails)}")
for f_ in fails:
    print("   FAIL", f_)


# ---------------------------------------------------------------- (1b) within-journal comparison (TMLR)
if len(sys.argv) > 4:
    TOUT = Path(sys.argv[4])
    T = json.loads((TOUT / "tmlr_results.json").read_text()); C = json.loads((TOUT / "tmlr_calibration.json").read_text())
    g14, g28, pt = T["groups"]["W14"], T["groups"]["W28"], T["passthrough"]
    acc = T["accounting"]; cs = T["common_support"]; sd = T["share_diffs"]; lh = T["late_hazard"]["pooled"]
    hz = {g: {x["day"]: x["hazard"] for x in T["hazard"][g]} for g in ("W14", "W28")}
    n0 = n
    hc = T["hazard_contrast"]; AE = T["action_editors"]; SEL = T["selection"]; fe = T["fe_log_chars"]["all"]["coef"]
    chk("tmlr papers 14", g14["n_papers"], 882, r"882 are same-day panels with a 14-day window (2,716 reviews)")
    chk("tmlr reviews 14", g14["n_reviews"], 2716, r"(2,716 reviews) and 233 with a 28-day window (710 reviews)")
    chk("tmlr papers 28", g28["n_papers"], 233, r"233 with a 28-day window (710 reviews)")
    chk("tmlr reviews 28", g28["n_reviews"], 710, r"233 with a 28-day window (710 reviews)")
    chk("tmlr panels", g14["n_papers"] + g28["n_papers"], 1115, r"In the 1,115 panels whose reviewers were all assigned on one day (3,426 reviews)")
    chk("tmlr reviews", g14["n_reviews"] + g28["n_reviews"], 3426, r"and on 3,426 reviews of a journal, TMLR,")
    chk("tmlr papers with due and reviews", acc["with_due_and_reviews"], 5560, r"Of 5,560 papers with a public due date and public reviews")
    chk("tmlr dropped panels", acc["dropped_reviews_before_under_review"]["papers"], 3, r"Three panels whose reviews are all dated before the recorded under-review date are left out")
    chk("tmlr all reviews", acc["reviews_all"], 17280, r"Over all 17,280 reviews of the journal")
    chk("tmlr lower bound late", 100 * T["full_sample_after_final_due"]["RL"]["share"], 38, r"at least 38\% are late", tol=0.5)
    chk("tmlr share same-day", (g14["n_papers"] + g28["n_papers"]) / acc["with_due_and_reviews"], 0.2, r"the sample is the fifth of the panels whose reviewers were all assigned on one day", tol=0.01)
    chk("tmlr median tau 14", g14["tau_q50"], 15.1, r"With a 14-day window the median review arrives 15.1 days after the assignment; with a 28-day window, 28.9 days after.")
    chk("tmlr median tau 28", g28["tau_q50"], 28.9, r"the median review arrives 15.1 days after the assignment when the window is 14 days and 28.9 days after when it is 28")
    chk("tmlr diff median", g28["tau_q50"] - g14["tau_q50"], 13.8, r"Fourteen more days of window go with a median delivery 13.8 days later")
    chk("tmlr pt q50", pt["q50"]["passthrough"], 0.99, r"a pass-through of 0.99, with a 95\% interval of 0.95 to 1.04 from resampling papers")
    chk("tmlr pt q50 lo", pt["q50"]["lo"], 0.95, r"a pass-through of 0.99; its 95\% interval, 0.95 to 1.04, measures sampling only")
    chk("tmlr pt q50 hi", pt["q50"]["hi"], 1.04, r"a pass-through of 0.99, with a 95\% interval of 0.95 to 1.04 from resampling papers")
    chk("tmlr abstract 13.8", g28["tau_q50"] - g14["tau_q50"], 13.8, r"At TMLR the median review arrives 13.8 days later under the longer window (0.99 per day)")
    chk("tmlr abstract 0.99", pt["q50"]["passthrough"], 0.99, r"13.8 days later under the longer window (0.99 per day)")
    chk("tmlr AE lo", AE["passthrough_q50_editor_bootstrap"]["lo"], 0.95, r"(0.95 to 1.05 from resampling action editors, and 1.00 within the 128 action editors who handled panels of both kinds")
    chk("tmlr AE hi", AE["passthrough_q50_editor_bootstrap"]["hi"], 1.05, r"(0.95 to 1.05 from resampling action editors")
    chk("tmlr AE within", AE["within_editor_weighted_median_diff_over_14"], 1.00, r"and 1.00 within the 128 action editors")
    chk("tmlr AE n", AE["n_editors_both_windows"], 128, r"within the 128 action editors who handled panels of both kinds")
    chk("tmlr s median 14", -g14["s_q50"], 0.6, r"the two medians fall 0.6 and 0.3 days \emph{after} it")
    chk("tmlr s median 28", -g28["s_q50"], 0.3, r"the two medians fall 0.6 and 0.3 days \emph{after} it")
    chk("tmlr elapsed 12-15 W14", 100 * T["elapsed_days_12_15"]["W14"]["share"], 35, r"35\% of the 14-day reviewers who have not yet delivered do so, against 4\% of the 28-day reviewers", tol=0.5)
    chk("tmlr elapsed 12-15 W28", 100 * T["elapsed_days_12_15"]["W28"]["share"], 4, r"against 4\% of the 28-day reviewers, who deliver at a similar rate (30\%) between days 26 and 29", tol=0.5)
    chk("tmlr elapsed 26-29 W28", 100 * T["elapsed_days_26_29"]["W28"]["share"], 30, r"(30\%) between days 26 and 29", tol=0.5)
    for g in ("W14", "W28"):
        chk(f"tmlr hazard day2 {g}", 100 * hz[g][2], 4.5, r"4--5\% two days before, 7\% the day before, and 19\% and 16\% in the final 24 hours", tol=0.5)
        chk(f"tmlr hazard day1 {g}", 100 * hz[g][1], 7, r"7\% the day before", tol=0.5)
        chk(f"tmlr hazard until a week before {g}", max(100 * v for d, v in hz[g].items() if d >= 7), 2, r"about 2\% per day or less until a week before it", tol=0.25)
    chk("tmlr hazard day0 W14", 100 * hz["W14"][0], 19, r"19\% and 16\% in the final 24 hours", tol=0.5)
    chk("tmlr hazard day0 W28", 100 * hz["W28"][0], 16, r"19\% and 16\% in the final 24 hours", tol=0.5)
    chk("tmlr pre-deadline hazard 14", 100 * hc["before_13_0"]["h14"], 3.4, r"Over the last fourteen days they average 3.4\% and 3.2\% per day, a difference of 0.3 points with an interval from $-0.7$ to 0.1")
    chk("tmlr pre-deadline hazard 28", 100 * hc["before_13_0"]["h28"], 3.2, r"they average 3.4\% and 3.2\% per day")
    chk("tmlr pre-deadline diff", -100 * hc["before_13_0"]["diff"], 0.3, r"a difference of 0.3 points with an interval from $-0.7$ to 0.1")
    chk("tmlr pre-deadline diff lo", 100 * hc["before_13_0"]["lo"], -0.7, r"with an interval from $-0.7$ to 0.1"); chk("tmlr pre-deadline diff hi", 100 * hc["before_13_0"]["hi"], 0.1, r"with an interval from $-0.7$ to 0.1")
    chk("tmlr first late day 14", 100 * hz["W14"][-1], 16.5, r"On the first late day 16.5\% of the pending 14-day reviews are delivered against 10.6\% of the 28-day ones")
    chk("tmlr first late day 28", 100 * hz["W28"][-1], 10.6, r"against 10.6\% of the 28-day ones")
    chk("tmlr late week 14", T["late_hazard"]["W14"]["1-7"]["hazard"], 0.13, r"over the first late week the daily hazard is 0.13 against 0.10 (a difference of 3.3 points; 1.5 to 5.0)")
    chk("tmlr late week 28", T["late_hazard"]["W28"]["1-7"]["hazard"], 0.10, r"the daily hazard is 0.13 against 0.10")
    chk("tmlr late week diff", -100 * hc["after_1_7"]["diff"], 3.3, r"(a difference of 3.3 points; 1.5 to 5.0)")
    chk("tmlr late week diff lo", -100 * hc["after_1_7"]["hi"], 1.5, r"(a difference of 3.3 points; 1.5 to 5.0)"); chk("tmlr late week diff hi", -100 * hc["after_1_7"]["lo"], 5.0, r"(a difference of 3.3 points; 1.5 to 5.0)")
    chk("tmlr a quarter more slowly", 1 - T["late_hazard"]["W28"]["1-7"]["hazard"] / T["late_hazard"]["W14"]["1-7"]["hazard"], 0.25, r"deliver about a quarter more slowly", tol=0.02)
    chk("tmlr pooled late hazards", lh["1-7"]["hazard"], 0.13, r"(0.13, 0.10 and 0.08 per day in the first week, the second, and the third and fourth, pooled)")
    chk("tmlr pt q25", pt["q25"]["passthrough"], 0.82, r"The lower quartile of the delivery date moves by 0.82 of the added days (0.65 to 0.92)")
    chk("tmlr pt q25 lo", pt["q25"]["lo"], 0.65, r"0.82 of the added days (0.65 to 0.92)")
    chk("tmlr pt q25 hi", pt["q25"]["hi"], 0.92, r"0.82 of the added days (0.65 to 0.92)")
    chk("tmlr pt q75", pt["q75"]["passthrough"], 1.12, r"the upper quartile moves by 1.12 (0.98 to 1.24)")
    chk("tmlr pt q75 lo", pt["q75"]["lo"], 0.98, r"the upper quartile moves by 1.12 (0.98 to 1.24)")
    chk("tmlr pt q75 hi", pt["q75"]["hi"], 1.24, r"the upper quartile moves by 1.12 (0.98 to 1.24)")
    chk("tmlr pt q90", pt["q90"]["passthrough"], 1.23, r"and the ninth decile by 1.23 (1.06 to 1.47)")
    chk("tmlr pt q90 lo", pt["q90"]["lo"], 1.06, r"the ninth decile by 1.23 (1.06 to 1.47)", tol=0.006)
    chk("tmlr pt q90 hi", pt["q90"]["hi"], 1.47, r"the ninth decile by 1.23 (1.06 to 1.47)")
    chk("tmlr first half", 100 * cs["share28_delivered_before_last14"], 12, r"Twelve per cent of the 28-day reviews are delivered in the first 14 days", tol=0.5)
    chk("tmlr KS p", cs["p_perm"], 0.0, r"a test on the days that the two windows share rejects the equality of the two distributions", tol=0.01)
    chk("tmlr early 28", 100 * sd["early_ge4"]["share28"], 26, r"more reviews delivered four days or more before the deadline (26\% against 18\%)", tol=0.5)
    chk("tmlr early 14", 100 * sd["early_ge4"]["share14"], 18, r"(26\% against 18\%)", tol=0.5)
    chk("tmlr final48 diff", 100 * sd["final48h"]["diff"], -4.1, r"fewer in the final 48 hours ($-4.1$ points; $-7.1$ to $-1.0$)", tol=0.05)
    chk("tmlr final48 diff lo", 100 * sd["final48h"]["lo"], -7.1, r"($-4.1$ points; $-7.1$ to $-1.0$)", tol=0.05); chk("tmlr final48 diff hi", 100 * sd["final48h"]["hi"], -1.0, r"($-4.1$ points; $-7.1$ to $-1.0$)", tol=0.05)
    chk("tmlr late diff", -100 * sd["late"]["diff"], 3.7, r"a late share lower by 3.7 points, the sign of P5, with an interval from $-8.1$ to 0.6", tol=0.05)
    chk("tmlr late diff lo", 100 * sd["late"]["lo"], -8.1, r"with an interval from $-8.1$ to 0.6", tol=0.05)
    chk("tmlr late diff hi", 100 * sd["late"]["hi"], 0.6, r"with an interval from $-8.1$ to 0.6", tol=0.05)
    chk("tmlr late>7 diff", 100 * sd["late_gt7"]["diff"], 3.6, r"and more reviews over a week late ($+3.6$ points; $-0.1$ to 7.6)", tol=0.05)
    chk("tmlr late>7 diff lo", 100 * sd["late_gt7"]["lo"], -0.1, r"($+3.6$ points; $-0.1$ to 7.6)", tol=0.05); chk("tmlr late>7 diff hi", 100 * sd["late_gt7"]["hi"], 7.6, r"($+3.6$ points; $-0.1$ to 7.6)", tol=0.05)
    chk("tmlr final48 14", 100 * g14["share_final48h"], 19, r"Only 19\% and 15\% of the reviews arrive in the final 48 hours, against 42--52\% at ICLR", tol=0.5)
    chk("tmlr final48 28", 100 * g28["share_final48h"], 15, r"Only 19\% and 15\% of the reviews arrive in the final 48 hours", tol=0.5)
    chk("tmlr late 14", 100 * g14["share_late"], 57, r"and 57\% and 53\% arrive after the deadline, against 16--37\%", tol=0.5)
    chk("tmlr late 28", 100 * g28["share_late"], 53, r"and 57\% and 53\% arrive after the deadline", tol=0.5)
    chk("tmlr late>7 of late 14", g14["share_late_gt7"] / g14["share_late"], 0.4, r"four in ten of the late reviews of the 14-day window", tol=0.02)
    chk("tmlr late>7 of late 28", g28["share_late_gt7"] / g28["share_late"], 0.5, r"and half of those of the 28-day window, are more than a week late", tol=0.02)
    chk("tmlr late hazard wk1", lh["1-7"]["hazard"], 0.13, r"After the deadline the daily hazard is 0.13 in the first week and 0.08--0.10 thereafter, against 0.29")
    chk("tmlr late hazard wk2", lh["8-14"]["hazard"], 0.10, r"0.08--0.10 thereafter")
    chk("tmlr late hazard wk3-4", lh["15-28"]["hazard"], 0.08, r"0.08--0.10 thereafter")
    chk("iclr late hazard", RES["E7"]["profile"]["mu"], 0.29, r"against 0.29 in the recent ICLR cohorts")
    lo_ = T["fe_log_chars_late_vs_on_time"]; f48 = T["fe_log_chars_final48h_vs_earlier_on_time"]; ce = T["fe_claims_supported_late_vs_on_time"]
    chk("tmlr fe final 24 h", fe["last24h"]["b"], 0.087, r"(final 24 hours $+0.087$ log points, standard error 0.047; late $-0.053$, 0.034)")
    chk("tmlr fe final 24 h se", fe["last24h"]["se"], 0.047, r"$+0.087$ log points, standard error 0.047")
    chk("tmlr fe late", fe["late"]["b"], -0.053, r"late $-0.053$, 0.034)"); chk("tmlr fe late se", fe["late"]["se"], 0.034, r"late $-0.053$, 0.034)")
    chk("tmlr fe none significant", max(abs(v["b"] / v["se"]) for v in fe.values()), 1.9, r"no timing bin differs significantly from early reviews", tol=0.06)
    chk("tmlr late length", lo_["b"], -0.098, r"Late reviews are 9\% shorter than on-time reviews taken together ($-0.098$; 0.021)")
    chk("tmlr late length se", lo_["se"], 0.021, r"($-0.098$; 0.021)")
    chk("tmlr late length pct", 100 * (1 - math.exp(lo_["b"])), 9, r"Late reviews are 9\% shorter than on-time reviews taken together", tol=0.5)
    chk("tmlr final48 length", f48["b"], 0.100, r"($+0.100$, 0.041, against earlier on-time reviews)")
    chk("tmlr final48 length se", f48["se"], 0.041, r"($+0.100$, 0.041, against earlier on-time reviews)")
    chk("tmlr claims", ce["b"], -0.007, r"($-0.007$; 0.018)")
    chk("tmlr claims se", ce["se"], 0.018, r"($-0.007$; 0.018)")
    chk("tmlr reviews 6% longer", 100 * (math.exp(T["log_chars_W28_minus_W14"]["diff"]) - 1), 6, r"under it more reviews arrive early and reviews are 6\% longer", tol=0.5)
    chk("tmlr same-day share regular", 100 * SEL["R"]["share_same_day"], 21, r"Same-day panels are 21\% of the regular submissions with reviews and 14\% of the long ones", tol=0.5)
    chk("tmlr same-day share long", 100 * SEL["L"]["share_same_day"], 14, r"21\% of the regular submissions with reviews and 14\% of the long ones", tol=0.5)
    chk("tmlr accept long same-day", 100 * SEL["L"]["same_day"]["accept_share_of_decided"], 64, r"long same-day papers were accepted less often than the other long papers (64\% against 74\%)", tol=0.5)
    chk("tmlr accept long other", 100 * SEL["L"]["other"]["accept_share_of_decided"], 74, r"(64\% against 74\%)", tol=0.5)
    chk("tmlr accept regular same-day", 100 * SEL["R"]["same_day"]["accept_share_of_decided"], 72.5, r"with no such difference among regular ones (72.5\% against 71.6\%; Table~S25)")
    chk("tmlr accept regular other", 100 * SEL["R"]["other"]["accept_share_of_decided"], 71.6, r"(72.5\% against 71.6\%; Table~S25)")
    chk("tmlr cal lo", C["range_journal_profile"][0], 25, r"gives journal optima of 25--49 days in the scenarios of Table~\ref{tab:tstar}, against 23--44")
    chk("tmlr cal hi", C["range_journal_profile"][1], 49, r"gives journal optima of 25--49 days")
    chk("cal conf lo", C["range_conference_profile"][0], 23, r"against 23--44: one to six days longer")
    chk("cal conf hi", C["range_conference_profile"][1], 44, r"against 23--44: one to six days longer")
    chk("tmlr cal shift lo", C["shift_range"][0], 1, r"one to six days longer")
    chk("tmlr cal shift hi", C["shift_range"][1], 6, r"one to six days longer")
    chk("tmlr cal shift0 lo", C["shift_range_zero_hazard_beyond_27"][0], -1, r"between one day shorter and two days longer if no report is delivered more than four weeks ahead")
    chk("tmlr cal shift0 hi", C["shift_range_zero_hazard_beyond_27"][1], 2, r"between one day shorter and two days longer")
    chk("tmlr mu", C["mu_tmlr"], 0.11, r"a late hazard of 0.11 instead of 0.29, and a shortfall for late reports only (0.094)")
    chk("tmlr Delta late", C["Delta_late_tmlr"], 0.094, r"a shortfall for late reports only (0.094)")
    chk("tmlr deadline-day hazard lo", 100 * hz["W28"][0], 16, r"a hazard of 16--19\% instead of 59\% on the deadline day", tol=0.5)
    chk("tmlr deadline-day hazard hi", 100 * hz["W14"][0], 19, r"a hazard of 16--19\% instead of 59\% on the deadline day", tol=0.5)
    chk("iclr deadline-day hazard", 100 * float(RES["E7"]["profile"]["h"]["0"]), 59, r"instead of 59\% on the deadline day", tol=0.5)
    chk("tmlr intro final48 range lo", 100 * g28["share_final48h"], 15, r"Only 15--19\% of the reviews arrive in the final 48 hours and 53--57\% arrive after the deadline", tol=0.5)
    chk("tmlr discussion median", g14["tau_q50"], 15, r"at TMLR a two-week window produces a median of 15 days and a ninth decile of 29", tol=0.5)
    chk("tmlr discussion q90", g14["tau_q90"], 29, r"a median of 15 days and a ninth decile of 29", tol=0.5)
    pooled48 = (g14["share_final48h"] * g14["n_reviews"] + g28["share_final48h"] * g28["n_reviews"]) / (g14["n_reviews"] + g28["n_reviews"])
    chk("tmlr conclusions final48 < 1/5", pooled48, 0.18, r"fewer than one review in five arrives in the final two days and more than half arrive late", tol=0.02)
    chk("tmlr conclusions slope", pt["q50"]["passthrough"], 0.99, r"a slope of 0.99 that describes where delivery is centred")
    chk("tmlr discussion 13.8", g28["tau_q50"] - g14["tau_q50"], 13.8, r"goes with a median delivery 13.8 days later and a late share lower by 3.7 points")
    chk("tmlr discussion shares", -100 * sd["final48h"]["diff"], 4.1, r"and 4.1 and 3.7 points for fourteen days at TMLR", tol=0.05)
    chk("dev_panel negatives", sum(1 for x in RES["E4"]["by_year"] if x["outcome"] == "dev_panel" and x["coef"] < 0), 32, r"the 36 coefficients are small and 32 of them negative")
    chk("dev_panel sig", sum(1 for x in RES["E4"]["by_year"] if x["outcome"] == "dev_panel" and abs(x["coef"] / x["se"]) > 1.96 and x["coef"] < 0), 3, r"the three that differ from zero at the 5\% level (2018 and 2020) are negative")
    h28 = {x["day"]: x for x in T["hazard"]["W28"]}
    b3 = [sum(h28[d]["events"] for d in (a, a + 1, a + 2)) / sum(h28[d]["at_risk"] for d in (a, a + 1, a + 2)) for a in (16, 19, 22, 25)]
    for w_, nm_ in ((plain, "text"), (letter, "letter")):
        frag = r"its daily hazard is 0.6--1.0\% in three-day averages, against the 0.2--0.6\% that the power law gives" if nm_ == "text" else r"0.6\% to 1.0\% per day in three-day averages, against the 0.2\% to 0.6\% that the power law gives"
        chk(f"{nm_} tmlr tail lo", 100 * min(b3), 0.6, frag, where=w_); chk(f"{nm_} tmlr tail hi", 100 * max(b3), 1.0, frag, where=w_)
        chk(f"{nm_} power law lo", 100 * 1.10 * 28 ** -1.85, 0.2, frag, where=w_); chk(f"{nm_} power law hi", 100 * 1.10 * 17 ** -1.85, 0.6, frag, where=w_)
    # letter
    chk("letter panels", g14["n_papers"] + g28["n_papers"], 1115, r"In the 1,115 panels whose reviewers were all assigned on one day (3,426 reviews)", where=letter)
    chk("letter pt", pt["q50"]["passthrough"], 0.99, r"and 28.9 days after under the 28-day window, a pass-through of 0.99", where=letter)
    chk("letter 13.8", g28["tau_q50"] - g14["tau_q50"], 13.8, r"Under the longer window the median review arrives 13.8 days later (0.99 days per day of window)", where=letter)
    chk("letter pt interval", pt["q50"]["lo"], 0.95, r"The interval of the pass-through, 0.95 to 1.04, therefore measures sampling only", where=letter)
    chk("letter elapsed", 100 * T["elapsed_days_12_15"]["W14"]["share"], 35, r"35\% of the pending 14-day reviewers deliver against 4\% of the 28-day reviewers", tol=0.5, where=letter)
    chk("letter length", f48["b"], 0.100, r"($+0.100$ log points, s.e.\ 0.041)", where=letter)
    chk("letter cal", C["range_journal_profile"][1], 49, r"the optima are 25--49 days against 23--44, which shows", where=letter)
    chk("letter shares", 100 * g14["share_late"], 57, r"15--19\% of the reviews arrive in the final 48 hours and 53--57\% late, against 42--52\% and 16--37\% at ICLR", tol=0.5, where=letter)
    print(f"(1b) TMLR claims checked: {n - n0}; failures in total: {len(fails)}")
    for f in fails:
        print("   FAIL", f)

# ---------------------------------------------------------------- (2) quotes
norm = lambda s: re.sub(r"\s+", " ", s.replace("\\rev{", "").replace("}", "").replace("{", "").replace("\\emph", "").replace("\\ldots", "")).strip()
P = norm(plain + " ".join(tf.read_text() for tf in (R1 / "tables").glob("*.tex")))
quotes = re.findall(r"``(.+?)''", letter or "", flags=re.S)
skip = ("ldots", "Timestamped Reports", "The main changes", "How to read this letter", "one journal experiment", "Final 48 hours", "fifteen-fold", "give the same results", "Your three requests", "How Long Should", "Corrections of our own", "0.085 units less valuable", "decision quality does not decline", "no significant increase in score deviation", "not less informative", "last 48 h", "One journal experiment", "the same review form", "The first use of OpenReview", "the first measurement of when")
bad = []
for q in quotes:
    if any(s in q for s in skip) or len(q.split()) < 3:
        continue
    if norm(q).rstrip(".,;:") not in P:
        bad.append(q[:90])
print(f"(2) letter quotes checked: {len(quotes)}; not found in manuscript: {len(bad)}" if letter is not None else "(2) response letter not present: its claims and quotes were skipped")
for b in bad:
    print("   MISSING:", b)

# ------------------------------------------------- (3) citations and labels
cited = set(k.strip() for g in re.findall(r"\\cite[a-z]*\{([^}]*)\}", tex) for k in g.split(","))
for tf in (R1 / "tables").glob("*.tex"):
    cited |= set(k.strip() for g in re.findall(r"\\cite[a-z]*\{([^}]*)\}", tf.read_text()) for k in g.split(","))
bib = set(re.findall(r"@\w+\{([^,]+),", (R1 / "references.bib").read_text()))
print(f"(3) cited {len(cited)}, bib {len(bib)}; cited not in bib: {sorted(cited - bib)}; in bib not cited: {sorted(bib - cited)}")
labels = set(re.findall(r"\\label\{([^}]*)\}", tex))
for tf in (R1 / "tables").glob("*.tex"):
    labels |= set(re.findall(r"\\label\{([^}]*)\}", tf.read_text()))
refs = set(re.findall(r"\\(?:eq)?ref\{([^}]*)\}", tex))
for tf in (R1 / "tables").glob("*.tex"):
    refs |= set(re.findall(r"\\(?:eq)?ref\{([^}]*)\}", tf.read_text()))
used_tables = set(re.findall(r"\\input\{tables/([^}]*)\}", tex))
print("    refs without label:", sorted(refs - labels), "| float labels never referenced:", sorted(l for l in labels - refs if l.startswith(("tab:", "fig:")) and not any(l == x for x in ("tab:contrasts", "tab:quality_year", "tab:sensitivity", "tab:tomography", "tab:validation"))))
print("    tables input:", sorted(used_tables))
