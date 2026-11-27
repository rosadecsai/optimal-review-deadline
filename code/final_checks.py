"""Checks added after the pre-submission review.

usage: python code/final_checks.py <repository root> <output folder> [reviews.csv(.gz)]

  F1  hazard invariance across ICLR cohorts: does the daily hazard at s days before the deadline vary with the window?
  F2  what the estimated hazards imply for the quantile slopes of the cohort table (Proposition 1 does not imply zero)
  F3  pass-through at ICLR on all reviews (late ones included) and without the 12-day cohort   [needs the review-level file]
  F4  compression elasticity from alternative contrasts, and the conference optimum as a function of the elasticity
  F5  the slopes of the cohort table with the published deadlines and windows                    [needs the review-level file]
  F6  score and deviation coefficients by rating scale (2018-2025 against 2026); reviews per timing bin   [idem]
  F7  the shortfall Delta(T) capped at the value of a report: effect on the optima
  F8  joint analysis with the functional form of the acceptance curve drawn as well
  F9  the late hazard of the reviewer model for three degrees of present bias (Proposition 2)
"""
import sys, json, math, os, itertools
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

REPO, OUT = Path(sys.argv[1]), Path(sys.argv[2])
REVIEWS = sys.argv[3] if len(sys.argv) > 3 else None
_argv = sys.argv; sys.argv = [_argv[0], str(REPO), str(OUT)]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import revision_analyses as RA            # definitions only (its analyses run under __main__)
sys.argv = _argv
RES = RA.RES
F = {}
e3 = pd.DataFrame(RES["E3"]["table"]).sort_values("year"); TW = {int(r.year): float(r["T"]) for _, r in e3.iterrows()}
years = sorted(TW); T = np.array([TW[y] for y in years])

def ols(x, y):
    n = len(x); sxx = ((x - x.mean()) ** 2).sum(); b = ((x - x.mean()) * (y - y.mean())).sum() / sxx
    a = y.mean() - b * x.mean(); se = math.sqrt(((y - a - b * x) ** 2).sum() / (n - 2) / sxx); tc = stats.t.ppf(0.975, n - 2)
    return dict(slope=float(b), se=float(se), lo=float(b - tc * se), hi=float(b + tc * se), n=int(n))

# ------------------------------------------------------------------ F1
h = pd.read_csv(REPO / "output" / "tables" / "hazard_by_day.csv")
h = h[h["year"].isin(years) & h["day_before"].between(0, 11)]
piv = h.pivot(index="day_before", columns="year", values="hazard").sort_index()
lg = np.log(piv)
dev = lg.sub(lg.mean(axis=1), axis=0)                       # log-hazard relative to the mean of the seven cohorts on that day
d_all = dev.mean(axis=0).reindex(years).to_numpy()
r_all = ols(T, d_all)
per_day = {int(s): ols(T, lg.loc[s].reindex(years).to_numpy()) for s in piv.index}
F["F1"] = dict(days="0-11 before the deadline (the days every cohort has)", cohorts=years, windows=[float(t) for t in T],
               pooled_log_hazard_slope=r_all, pooled_pct_per_day=100 * (math.exp(r_all["slope"]) - 1),
               pooled_pct_lo=100 * (math.exp(r_all["lo"]) - 1), pooled_pct_hi=100 * (math.exp(r_all["hi"]) - 1),
               near=ols(T, dev.loc[0:3].mean(axis=0).reindex(years).to_numpy()), far=ols(T, dev.loc[4:11].mean(axis=0).reindex(years).to_numpy()),
               n_days_slope_significant=int(sum(1 for v in per_day.values() if v["lo"] > 0 or v["hi"] < 0)), per_day=per_day,
               hazard_range_by_day={int(s): [float(piv.loc[s].min()), float(piv.loc[s].max())] for s in piv.index})
print("F1 pooled slope of log hazard on T:", {k: round(v, 4) for k, v in r_all.items()}, "| % per day:", round(F["F1"]["pooled_pct_per_day"], 2),
      "| days with a significant slope:", F["F1"]["n_days_slope_significant"], "of", len(per_day))

# ------------------------------------------------------------------ F2
def on_time_quantiles(Tw, hz, ps=(0.25, 0.5, 0.75, 0.9)):
    """quantiles of days before the deadline among on-time reviews (uniform within the day), window Tw, hazards hz(s)"""
    surv = 1.0; mass = {}
    for s in range(int(Tw) - 1, -1, -1):
        hh = hz(s); mass[s] = surv * hh; surv *= 1 - hh
    tot = sum(mass.values()); out = []
    for p in ps:
        acc = 0.0
        for s in range(0, int(Tw)):
            if acc + mass[s] / tot >= p:
                out.append(s + (p - acc) / (mass[s] / tot)); break
            acc += mass[s] / tot
    return out
obs = {r["stat"]: r for r in RA.json.loads((OUT / "revision_results.json").read_text())["A1"]["table"]}
F["F2"] = {}
for tail in ("power", "exponential", "floor"):
    Q = np.array([on_time_quantiles(t, lambda s: RA.hazard(s, tail)) for t in T])
    F["F2"][tail] = {f"q{int(100*p)}": dict(implied_slope=ols(T, Q[:, i])["slope"], observed_slope=obs[f"s_q{int(100*p)}"]["slope"],
                                              observed_lo=obs[f"s_q{int(100*p)}"]["lo"], observed_hi=obs[f"s_q{int(100*p)}"]["hi"],
                                              implied_at_12=float(Q[list(T).index(12.0), i]), implied_at_21=float(Q[list(T).index(21.0), i]))
                     for i, p in enumerate((0.25, 0.5, 0.75, 0.9))}
for tail in ("power", "exponential", "floor"):
    atom_, late_, _, _ = RA.PRE[tail]; ix = [int(np.where(RA.TS == t)[0][0]) for t in T]
    F["F2"][tail]["share_final48h"] = dict(implied_slope=ols(T, 100 * atom_[ix])["slope"], observed_slope=100 * obs["share_atom"]["slope"])
    F["F2"][tail]["share_late"] = dict(implied_slope=ols(T, 100 * late_[ix])["slope"], observed_slope=100 * obs["share_late"]["slope"])
print("F2 implied vs observed slopes (power tail):", {k: (round(v["implied_slope"], 3), round(v["observed_slope"], 3)) for k, v in F["F2"]["power"].items()})

# ------------------------------------------------------------------ F3
on = ols(T, e3["s_q50"].to_numpy(float))
F["F3"] = dict(on_time=dict(passthrough=1 - on["slope"], lo=1 - on["hi"], hi=1 - on["lo"]))
k = [i for i, y in enumerate(years) if y != 2023]
o6 = ols(T[k], e3["s_q50"].to_numpy(float)[k])
F["F3"]["on_time_without_2023"] = dict(passthrough=1 - o6["slope"], lo=1 - o6["hi"], hi=1 - o6["lo"], n=6)
if REVIEWS:
    df = pd.read_csv(REVIEWS, low_memory=False, usecols=["year", "t_before", "in_sample"])
    df = df[df["in_sample"].astype(bool) & df["year"].isin(years)]
    g = df.groupby("year")["t_before"]
    med_on = df[df["t_before"] >= 0].groupby("year")["t_before"].median().reindex(years).to_numpy()
    assert np.allclose(med_on, e3["s_q50"].to_numpy(float), atol=1e-6), "on-time medians differ from the archived ones"
    res = {}
    for p in (0.25, 0.5, 0.75):
        qa = g.quantile(p).reindex(years).to_numpy(); a = ols(T, qa); a6 = ols(T[k], qa[k])
        sxx = ((T - T.mean()) ** 2).sum(); bp = np.array([((T - T.mean()) * (qa[list(pp)] - qa.mean())).sum() / sxx for pp in itertools.permutations(range(len(T)))])
        res[f"q{int(100*p)}"] = dict(passthrough=1 - a["slope"], lo=1 - a["hi"], hi=1 - a["lo"], slope_days_before=a["slope"], p_perm=float((np.abs(bp) >= abs(a["slope"]) - 1e-15).mean()),
                                     by_cohort={int(y): float(v) for y, v in zip(years, qa)}, without_2023=dict(passthrough=1 - a6["slope"], lo=1 - a6["hi"], hi=1 - a6["lo"]))
    F["F3"]["all_reviews"] = res
    print("F3 pass-through of the median: on-time", {k_: round(v, 3) for k_, v in F["F3"]["on_time"].items()}, "| all reviews", {k_: round(v, 3) for k_, v in res["q50"].items() if k_ in ("passthrough", "lo", "hi")},
          "| on-time without 2023", {k_: round(v, 3) for k_, v in F["F3"]["on_time_without_2023"].items()})

# ------------------------------------------------------------------ F4
by = pd.DataFrame(RES["E4"]["by_year"]); lw = by[(by["outcome"] == "log_words") & (by["bin"] == "last24h")].set_index("year")
def pool(ys):
    x = lw.loc[list(ys)]; w = 1 / x["se"] ** 2
    return float((w * x["coef"]).sum() / w.sum()), float(math.sqrt(1 / w.sum()))
g23 = float(lw.loc[2023, "coef"]); contrasts = {}
for name, ys in (("2024-2026 (as in the calibration)", (2024, 2025, 2026)), ("2020-2022 (before the 2024 form)", (2020, 2021, 2022)), ("2018-2022", (2018, 2019, 2020, 2021, 2022)), ("2022 (adjacent, same scores)", (2022,))):
    est, se = pool(ys); Tref = float(np.mean([TW[y] for y in ys if y in TW]))
    eta = math.log(g23 / est) / math.log(Tref / 12.0) if est < 0 else float("nan")
    contrasts[name] = dict(gap=est, se=se, T_ref=Tref, eta=eta)
single = {int(y): dict(gap=float(lw.loc[y, "coef"]), se=float(lw.loc[y, "se"]), T=TW[y], eta=(math.log(g23 / float(lw.loc[y, "coef"])) / math.log(TW[y] / 12.0))) for y in years if y != 2023}
F["F4"] = dict(gap_2023=g23, se_2023=float(lw.loc[2023, "se"]), contrasts=contrasts, single_cohorts=single, Delta_final48h=RA.DA0, Delta_late=RA.DL0)
print("F4 eta against single cohorts:", {y: round(v["eta"], 2) for y, v in single.items()})
etas = [0.0, 0.25, 0.5, 0.7, 1.0, 1.3, 1.6, 2.0]; rows = []
for eta in etas:
    r = dict(eta=eta)
    for delta in (0.0025, 0.005, 0.01):
        V = RA.value(delta, 0.0, RA.DA0, RA.DL0, eta, np.ones((1, len(RA.TS))))
        r[str(delta)] = int(RA.TS[np.argmax(V[0])])
    rows.append(r)
F["F4"]["conference_optimum_by_eta"] = rows
print("F4 eta by contrast:", {k_: round(v["eta"], 2) for k_, v in contrasts.items()}); print(pd.DataFrame(rows).to_string(index=False))

# ------------------------------------------------------------------ F5-F6 (review-level file)
PUB = {2020: 20, 2021: 16, 2022: 18, 2023: 12, 2024: 21, 2025: 20, 2026: 22}     # windows of the reviewer guides
if REVIEWS:
    df = pd.read_csv(REVIEWS, low_memory=False, usecols=["year", "forum", "t_before", "t_before_official", "in_sample", "bin", "score", "dev_panel", "n_reviews"])
    df = df[df["in_sample"].astype(bool)]
    Tp = np.array([PUB[y] for y in years], float); res = {}
    for col, Tw, lab in (("t_before", T, "effective"), ("t_before_official", Tp, "published")):
        d = df[df["year"].isin(years)]; on = d[d[col] >= 0]; g = on.groupby("year")[col]; r = {}
        for q in (0.25, 0.5, 0.75, 0.9):
            y = g.quantile(q).reindex(years).to_numpy(); o = ols(Tw, y)
            r[f"q{int(100*q)}"] = dict(o, by_cohort=[float(v) for v in y])
        r["passthrough_median"] = dict(est=1 - r["q50"]["slope"], lo=1 - r["q50"]["hi"], hi=1 - r["q50"]["lo"])
        for nm, fn in (("share_final48h", lambda x: ((x >= 0) & (x < 2)).mean()), ("share_late", lambda x: (x < 0).mean())):
            y = d.groupby("year")[col].apply(fn).reindex(years).to_numpy(); r[nm] = dict(ols(Tw, 100 * y), by_cohort=[float(100 * v) for v in y])
        res[lab] = r
    assert abs(res["effective"]["q50"]["slope"] - obs["s_q50"]["slope"]) < 1e-9, "effective-deadline slope differs from the archived one"
    F["F5"] = dict(published_windows=[PUB[y] for y in years], **res)
    print("F5 slopes with the published deadlines and windows:", {k_: (round(v["slope"], 3), round(v["lo"], 3), round(v["hi"], 3)) for k_, v in res["published"].items() if "slope" in v},
          "| pass-through of the median", {k_: round(v, 3) for k_, v in res["published"]["passthrough_median"].items()})

    def fe(d, yname):
        """manuscript fixed effects, four timing dummies (reference: early), standard errors clustered by manuscript"""
        d = d.dropna(subset=[yname]); d = d[d.groupby("forum")[yname].transform("size") >= 2]
        bins = ["mid", "last48h", "last24h", "late"]
        X = np.column_stack([(d["bin"] == b).to_numpy(float) for b in bins]); y = d[yname].to_numpy(float)
        codes = pd.factorize(d["forum"])[0]; G = codes.max() + 1; cnt = np.bincount(codes, minlength=G)
        y = y - (np.bincount(codes, weights=y, minlength=G) / cnt)[codes]
        for j in range(X.shape[1]):
            X[:, j] -= (np.bincount(codes, weights=X[:, j], minlength=G) / cnt)[codes]
        XtXi = np.linalg.inv(X.T @ X); b = XtXi @ X.T @ y; e = y - X @ b
        S = np.column_stack([np.bincount(codes, weights=X[:, j] * e, minlength=G) for j in range(X.shape[1])])
        n, k = X.shape; V = XtXi @ (S.T @ S) @ XtXi * (G / (G - 1)) * ((n - 1) / (n - k - G))
        return dict(n=int(n), manuscripts=int(G), coef={bb: dict(b=float(b[i]), se=float(np.sqrt(V[i, i]))) for i, bb in enumerate(bins)},
                    bin_n={bb: int((d["bin"] == bb).sum()) for bb in ["early"] + bins})
    f6 = {}
    for yname in ("score", "dev_panel"):
        f6[yname] = {"2018-2026": fe(df, yname), "2018-2025 (1-10 scale)": fe(df[df["year"] <= 2025], yname), "2026 (0-10 scale in steps of two)": fe(df[df["year"] == 2026], yname)}
    pooled = pd.DataFrame(RES["E4"]["pooled"]); chk = pooled[(pooled["outcome"] == "score") & (pooled["bin"] == "last24h")].iloc[0]
    assert abs(f6["score"]["2018-2026"]["coef"]["last24h"]["b"] - float(chk["coef"])) < 2e-3, "pooled score coefficient not reproduced"
    F["F6"] = f6
    print("F6 final-24-hour score coefficient:", {k_: (round(v["coef"]["last24h"]["b"], 3), round(v["coef"]["last24h"]["se"], 3)) for k_, v in f6["score"].items()},
          "| deviation:", {k_: (round(v["coef"]["last24h"]["b"], 3), round(v["coef"]["last24h"]["se"], 3)) for k_, v in f6["dev_panel"].items()})
    print("   reviews by timing bin (regression sample):", f6["score"]["2018-2026"]["bin_n"], "n", f6["score"]["2018-2026"]["n"])

# ------------------------------------------------------------------ F7: the shortfall capped at the value of a report
def value_capped(delta, L, Da, Dl, eta, acc, tail="power", mu=RA.MU0, rho=1.0, T0=28, cap=1.0):
    delta, L, Da, Dl, eta, mu, rho = [np.atleast_1d(np.asarray(x, dtype=float)) for x in (delta, L, Da, Dl, eta, mu, rho)]
    atom, late, _, _ = RA.PRE[tail]; comp = (RA.T_REF / RA.TS[None, :]) ** eta[:, None]
    Eq = 1.0 - np.minimum(cap, Da[:, None] * comp) * atom[None, :] - np.minimum(cap, Dl[:, None] * comp) * late[None, :]
    d = RA.delay(tail, mu); i0 = int(np.where(RA.TS == T0)[0][0]); d = d[:, [i0]] + rho[:, None] * (d - d[:, [i0]])
    return acc * (Eq - delta[:, None] * d) - (1 - acc) * L[:, None]
dr = np.load(OUT / "A7_draws.npz"); n = len(dr["delta"]); Tc0 = np.zeros(n, int); Tj0 = np.zeros(n, int); Tc1 = np.zeros(n, int); Tj1 = np.zeros(n, int)
for ti, tname in enumerate(RA.TAILS):
    sel = dr["tail"] == ti; one = np.ones((sel.sum(), len(RA.TS))); aj = RA.acc_exp(dr["lam"][sel])
    a = (dr["delta"][sel], 0.0, dr["Delta"][sel], dr["Delta"][sel], dr["eta"][sel], one); b = (dr["delta"][sel], dr["L"][sel], dr["Delta"][sel], dr["Delta"][sel], dr["eta"][sel], aj)
    Tc0[sel] = RA.tstar(RA.value(*a, tail=tname, mu=dr["mu"][sel])); Tj0[sel] = RA.tstar(RA.value(*b, tail=tname, mu=dr["mu"][sel], rho=dr["rho"][sel]))
    Tc1[sel] = RA.tstar(value_capped(*a, tail=tname, mu=dr["mu"][sel])); Tj1[sel] = RA.tstar(value_capped(*b, tail=tname, mu=dr["mu"][sel], rho=dr["rho"][sel]))
assert (Tc0 == dr["Tc"]).all() and (Tj0 == dr["Tj"]).all(), "joint draws not reproduced"
qs = (5, 25, 50, 75, 95)
tab10 = []
for eta in (0.0, 1.3):
    for delta in (0.0025, 0.005, 0.01):
        one = np.ones((1, len(RA.TS))); sel60 = RA.TS <= 60
        t0 = int(RA.TS[sel60][np.argmax(RA.value(delta, 0.0, RA.DA0, RA.DL0, eta, one)[0][sel60])]); t1 = int(RA.TS[sel60][np.argmax(value_capped(delta, 0.0, RA.DA0, RA.DL0, eta, one)[0][sel60])])
        tab10.append(dict(eta=eta, delta=delta, uncapped=t0, capped=t1))
dl = lambda T, D, eta: D * (RA.T_REF / T) ** eta
F["F7"] = dict(cap=1.0, conference=dict(q_uncapped={q: float(np.percentile(Tc0, q)) for q in qs}, q_capped={q: float(np.percentile(Tc1, q)) for q in qs}, share_changed=float((Tc0 != Tc1).mean()), max_change=int(np.abs(Tc0 - Tc1).max())),
               journal=dict(q_uncapped={q: float(np.percentile(Tj0, q)) for q in qs}, q_capped={q: float(np.percentile(Tj1, q)) for q in qs}, share_changed=float((Tj0 != Tj1).mean()), max_change=int(np.abs(Tj0 - Tj1).max())),
               share_draws_cap_binds_at_T3=float((dl(3.0, dr["Delta"], dr["eta"]) > 1).mean()), share_draws_cap_binds_at_T7=float((dl(7.0, dr["Delta"], dr["eta"]) > 1).mean()),
               share_draws_cap_binds_at_T12=float((dl(12.0, dr["Delta"], dr["eta"]) > 1).mean()),
               table10_conference=tab10, Delta_T={str(T): dl(T, RA.DA0, 1.3) for T in (3, 5, 7, 12, 21)}, Delta_T_upper_corner={str(T): dl(T, 0.30, 2.0) for T in (3, 7, 11, 12)},
               reference_window="20.5 = mean window of the 2024-2026 cohorts (21, 20, 21 days) on which the profile and the reference shortfall are estimated")
print("F7 shortfall capped at 1: conference", F["F7"]["conference"], "| journal", F["F7"]["journal"], "| Table 10 conference", tab10)

# ------------------------------------------------------------------ F8: joint analysis with the functional form of the acceptance curve drawn too
a28, a45 = 0.641, 0.676; rr_ = a45 / a28; kk = 28 * 45 * (rr_ - 1) / (45 - 28 * rr_); ab_h = a28 * (28 + kk) / 28
b_log = (a45 - a28) / math.log(45 / 28); c_log = a28 - b_log * math.log(28); b_lin = (a45 - a28) / (45 - 28); c_lin = a28 - b_lin * 28
TSf = RA.TS.astype(float)
FORMS = {"exponential": RA.ABAR0 * (1 - np.exp(-RA.LAM0 * TSf)), "hyperbolic": ab_h * TSf / (TSf + kk), "logarithmic": np.clip(c_log + b_log * np.log(TSf), 0, 1), "linear": np.clip(c_lin + b_lin * TSf, 0, 1)}
g = np.random.Generator(np.random.PCG64(RA.SEED + 8)); form = g.integers(0, 4, n); Tf = np.zeros(n, int); fn = list(FORMS)
for ti, tname in enumerate(RA.TAILS):
    for fi, fname in enumerate(fn):
        sel = (dr["tail"] == ti) & (form == fi)
        if sel.any():
            acc = np.repeat(FORMS[fname][None, :], sel.sum(), axis=0)
            Tf[sel] = RA.tstar(RA.value(dr["delta"][sel], dr["L"][sel], dr["Delta"][sel], dr["Delta"][sel], dr["eta"][sel], acc, tail=tname, mu=dr["mu"][sel], rho=dr["rho"][sel]))
gm = Tf.mean(); between = sum((form == fi).mean() * (Tf[form == fi].mean() - gm) ** 2 for fi in range(4))
F["F8"] = dict(n=int(n), q={q: float(np.percentile(Tf, q)) for q in qs}, share_at_90=float((Tf == 90).mean()), share_le7=float((Tf <= 7).mean()), share_lt21=float((Tf < 21).mean()), share_gt60=float((Tf > 60).mean()),
               by_form={fname: dict(q={q: float(np.percentile(Tf[form == fi], q)) for q in qs}, share_at_90=float((Tf[form == fi] == 90).mean()), share_le7=float((Tf[form == fi] <= 7).mean())) for fi, fname in enumerate(fn)},
               variance_share_form=float(between / Tf.var()), note="each of the four curves through the two arms with probability 1/4, in place of the curvature of the exponential; other inputs as in the joint draws")
print("F8 journal optimum with the functional form drawn:", F["F8"]["q"], "at 90:", round(F["F8"]["share_at_90"], 3), "<=7:", round(F["F8"]["share_le7"], 3), "| variance share of the form:", round(F["F8"]["variance_share_form"], 3))
print("   by form:", {k_: v["q"] for k_, v in F["F8"]["by_form"].items()})

# ------------------------------------------------------------------ F9: the late hazard in the reviewer model, by present bias
sys.path.insert(0, str(REPO / "code"))
from model import ReviewerModel
F["F9"] = dict(setting="log-normal costs of unit mean (sigma 0.8), lateness cost 0.25 per day, last admissible day 120 days after the deadline", late_hazard={})
for beta in (1.0, 0.9, 0.7):
    hh = ReviewerModel(sigma=0.8, beta=beta, late_pen=0.25, late_growth=0.0, max_late=120).hazard()
    F["F9"]["late_hazard"][str(beta)] = {str(d): float(hh[-d]) for d in (1, 5, 10, 20, 40)}
print("F9 late hazard by days late (1, 5, 10, 20, 40):", {b: [round(v, 3) for v in x.values()] for b, x in F["F9"]["late_hazard"].items()})

(OUT / "final_checks.json").write_text(json.dumps(F, indent=1, default=float))
print("written", OUT / "final_checks.json")
