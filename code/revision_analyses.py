"""Revision analyses (R1) for 'How Long Should Reviewers Get?'.

Everything here is computed from the archived outputs of the original pipeline
(output/results.json and data/processed/papers.csv.gz of the public repository);
no review-level record is needed. Seed: PCG64(20261126).

  A1  cohort-level inference for the Parkinson slopes (seven cohorts)
  A2  reproduction of the submitted optimum (check of this re-implementation)
  A3  acceptance curve: alternative functional forms through the two RCT arms
  A4  partial anchoring (pass-through rho) and the generalised closed form
  A5  tail of the hazard profile: alternative extrapolations
  A6  shortfall Delta from 0 to 0.30
  A7  joint (global) sensitivity of T*: Latin hypercube, first-order indices
  A8  shape of V(T): unimodality across the joint draws
  A9  length gap before/after the release of ChatGPT; equivalence test
  A10 agreement by timing pair: pooled early-early vs late-late contrast
  A11 number of reviews per manuscript (completion proxy)
  A12 closed form at the calibrated curvature and other numbers quoted in the text

Usage:  python code/revision_analyses.py <repository root> <output folder>
"""
from __future__ import annotations

import gzip
import io
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import brentq
from scipy.special import lambertw

REPO = Path(sys.argv[1])
OUT = Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
RES = json.loads((REPO / "output" / "results.json").read_text())
SEED = 20261126
OUTJ: dict = {}

# ---------------------------------------------------------------- profile ---
PROF = RES["E7"]["profile"]
H_EMP = {int(k): float(v) for k, v in PROF["h"].items()}
MU0, PL_A, PL_B = float(PROF["mu"]), float(PROF["pl_a"]), float(PROF["pl_b"])
DA0, DL0 = float(RES["E7"]["Delta_atom"]), float(RES["E7"]["Delta_late"])
S_EXT, MAX_LATE, T_REF = 15, 60, 20.5
TS = np.arange(3, 91)

# exponential tail fitted on the same range as the power law (s = 2..15)
_ss = np.arange(2, S_EXT + 1)
_hh = np.array([H_EMP[int(k)] for k in _ss])
EXP_K, EXP_C = np.polyfit(_ss, np.log(_hh), 1)
EXP_K = -EXP_K


def hazard(s: int, tail: str) -> float:
    """Pre-deadline daily hazard, s = days before the deadline (s >= 0)."""
    if s < S_EXT:
        return H_EMP[s]
    pl = min(0.999, PL_A * (1 + s) ** (-PL_B))
    if tail == "power":
        return pl
    if tail == "exponential":
        return float(np.exp(EXP_C - EXP_K * s))
    if tail == "truncated":       # nobody writes more than 30 days ahead
        return pl if s <= 30 else 0.0
    if tail == "floor":           # hazard stops decaying at the edge of the data
        return min(0.999, PL_A * (1 + S_EXT) ** (-PL_B))
    raise ValueError(tail)


TAILS = ["power", "exponential", "truncated", "floor"]


def pre_moments(tail: str, ts=TS):
    """For each T: share in the final 48 h, mass reaching the deadline unwritten
    (late share), sum of s*p over on-time days, and share reaching the final 48 h."""
    atom, late, spre, reach = [], [], [], []
    for T in ts:
        surv, a, sp, r = 1.0, 0.0, 0.0, None
        for s in range(int(T) - 1, -1, -1):
            if s == 1:
                r = surv
            h = hazard(s, tail)
            p = surv * h
            if s < 2:
                a += p
            sp += s * p
            surv *= 1 - h
        atom.append(a); late.append(surv); spre.append(sp); reach.append(r if r is not None else 1.0)
    return np.array(atom), np.array(late), np.array(spre), np.array(reach)


def mean_days_late(mu):
    """Mean days late (d = 1..60) of a geometric tail with daily hazard mu and
    forced delivery on day 60, as in the submitted DeliveryProfile."""
    mu = np.asarray(mu, dtype=float)
    d = np.arange(1, MAX_LATE)
    core = (d[None, :] * mu[..., None] * (1 - mu[..., None]) ** (d[None, :] - 1)).sum(-1)
    return core + MAX_LATE * (1 - mu) ** (MAX_LATE - 1)


PRE = {t: pre_moments(t) for t in TAILS}


def delay(tail: str, mu):
    """Expected delivery date (days after assignment) for every T in TS."""
    atom, late, spre, _ = PRE[tail]
    mu = np.atleast_1d(np.asarray(mu, dtype=float))
    return TS[None, :] - (spre[None, :] - late[None, :] * mean_days_late(mu)[:, None] + 0.5)


def value(delta, L, Da, Dl, eta, acc, tail="power", mu=MU0, rho=1.0, T0=28):
    """Editor's value V(T) on TS for vectors of parameters (one row per draw).
    acc: array (n, len(TS)) of acceptance probabilities (ones for a conference)."""
    delta, L, Da, Dl, eta, mu, rho = [np.atleast_1d(np.asarray(x, dtype=float)) for x in (delta, L, Da, Dl, eta, mu, rho)]
    atom, late, _, _ = PRE[tail]
    comp = (T_REF / TS[None, :]) ** eta[:, None]
    Eq = 1.0 - Da[:, None] * comp * atom[None, :] - Dl[:, None] * comp * late[None, :]
    d = delay(tail, mu)
    i0 = int(np.where(TS == T0)[0][0])
    d = d[:, [i0]] + rho[:, None] * (d - d[:, [i0]])
    return acc * (Eq - delta[:, None] * d) - (1 - acc) * L[:, None]


def acc_exp(lam, a28=0.641):
    lam = np.atleast_1d(np.asarray(lam, dtype=float))
    abar = a28 / (1 - np.exp(-28 * lam))
    return abar[:, None] * (1 - np.exp(-lam[:, None] * TS[None, :]))


def calibrate_exp(a28=0.641, a45=0.676):
    f = lambda lam: (1 - math.exp(-45 * lam)) / (1 - math.exp(-28 * lam)) - a45 / a28
    lam = brentq(f, 1e-4, 1.0)
    return a28 / (1 - math.exp(-28 * lam)), lam


ABAR0, LAM0 = calibrate_exp()


def tstar(V):
    return TS[np.argmax(V, axis=1)]


# ============================================================ A1 cohorts ====
def a1():
    tab = pd.DataFrame(RES["E3"]["table"]).sort_values("year")
    T = tab["T"].to_numpy(float)
    rows = []
    perms = list(itertools.permutations(range(len(T))))
    for stat, label in [("s_q25", "Q25 days before deadline"), ("s_q50", "Median days before deadline"),
                        ("s_q75", "Q75 days before deadline"), ("s_q90", "Q90 days before deadline"),
                        ("share_atom", "Share in final 48 h"), ("share_late", "Share late"), ("gamma", "Anchoring rate")]:
        y = tab[stat].to_numpy(float)
        n = len(T)
        sxx = ((T - T.mean()) ** 2).sum()
        b = ((T - T.mean()) * (y - y.mean())).sum() / sxx
        a = y.mean() - b * T.mean()
        res = y - a - b * T
        se = math.sqrt((res ** 2).sum() / (n - 2) / sxx)
        tcrit = stats.t.ppf(0.975, n - 2)
        # exact permutation test of slope = 0 (all 7! orderings of the cohort statistic)
        bp = np.array([((T - T.mean()) * (y[list(p)] - y.mean())).sum() / sxx for p in perms])
        p_perm = float((np.abs(bp) >= abs(b) - 1e-15).mean())
        loo = [np.polyfit(np.delete(T, i), np.delete(y, i), 1)[0] for i in range(n)]
        rows.append(dict(stat=stat, label=label, slope=b, se=se, lo=b - tcrit * se, hi=b + tcrit * se, t=b / se,
                         p_t=float(2 * stats.t.sf(abs(b / se), n - 2)), p_perm=p_perm, loo_min=min(loo), loo_max=max(loo),
                         loo_drop_min=int(tab["year"].iloc[int(np.argmin(loo))]), loo_drop_max=int(tab["year"].iloc[int(np.argmax(loo))])))
    df = pd.DataFrame(rows)
    # leverage of the 12-day cohort; correlation of T with calendar year
    lev = (1 / len(T) + (T - T.mean()) ** 2 / ((T - T.mean()) ** 2).sum())
    extra = dict(n=len(T), tcrit=float(stats.t.ppf(0.975, len(T) - 2)), leverage={int(y): float(l) for y, l in zip(tab["year"], lev)},
                 corr_T_year=float(np.corrcoef(T, tab["year"])[0, 1]), T=[float(x) for x in T], years=[int(x) for x in tab["year"]])
    print("\n== A1 cohort-level inference (n=7, t with 5 df) ==")
    print(df.round(4).to_string(index=False)); print(extra)
    OUTJ["A1"] = dict(table=df.to_dict(orient="records"), **extra)
    df.to_csv(OUT / "A1_cohort_slopes.csv", index=False)


# ======================================================= A2 reproduction ====
def a2():
    rows = []
    ok = True
    for sc in RES["E7"]["tstar"]:
        if sc["conference"]:
            acc = np.ones((1, len(TS)))
        else:
            acc = ABAR0 * (1 - np.exp(-LAM0 * TS[None, :]))
        V = value(sc["delta"], sc["L"], DA0, DL0, sc["eta"], acc)
        sel = TS <= 60
        mine = int(TS[sel][np.argmax(V[0][sel])])
        ok &= mine == sc["T_star"]
        rows.append(dict(scenario=sc["scenario"], submitted=sc["T_star"], reproduced=mine, V_sub=sc["V_star"], V_rep=float(V[0][sel].max())))
    df = pd.DataFrame(rows)
    print("\n== A2 reproduction of submitted Table 7 ==")
    print(df.round(4).to_string(index=False)); print("all identical:", ok, "| abar, lam:", round(ABAR0, 4), round(LAM0, 4))
    OUTJ["A2"] = dict(all_identical=bool(ok), max_abs_V_diff=float((df["V_sub"] - df["V_rep"]).abs().max()), abar=ABAR0, lam_a=LAM0)


# ===================================================== A3 acceptance forms ==
def a3():
    a28, a45 = 0.641, 0.676
    # hyperbolic a = abar*T/(T+k)
    r = a45 / a28
    k = 28 * 45 * (r - 1) / (45 - 28 * r)
    ab_h = a28 * (28 + k) / 28
    # logarithmic a = c + b ln T
    b_log = (a45 - a28) / math.log(45 / 28); c_log = a28 - b_log * math.log(28)
    # linear (chord through the two arms): upper envelope of every concave curve outside [28,45]
    b_lin = (a45 - a28) / (45 - 28); c_lin = a28 - b_lin * 28
    forms = {
        "exponential (submitted)": ABAR0 * (1 - np.exp(-LAM0 * TS)),
        "hyperbolic": ab_h * TS / (TS + k),
        "logarithmic": np.clip(c_log + b_log * np.log(TS), 0, 1),
        "linear (chord)": np.clip(c_lin + b_lin * TS, 0, 1),
    }
    params = dict(exponential=dict(abar=ABAR0, lam=LAM0), hyperbolic=dict(abar=ab_h, k=k), logarithmic=dict(c=c_log, b=b_log), linear=dict(c=c_lin, b=b_lin))
    rows = []
    for name, a in forms.items():
        for delta in (0.0025, 0.005, 0.01):
            for eta in (0.0, 1.3):
                V = value(delta, 0.3, DA0, DL0, eta, a[None, :])[0]
                i = int(np.argmax(V)); Ts = int(TS[i])
                rows.append(dict(form=name, delta=delta, eta=eta, T_star=Ts, a7=float(a[TS == 7][0]), a14=float(a[TS == 14][0]), a21=float(a[TS == 21][0]),
                                 loss14=float(100 * (V[i] - V[TS == 14][0])), loss42=float(100 * (V[i] - V[TS == 42][0]))))
    df = pd.DataFrame(rows)
    print("\n== A3 acceptance forms through a(28)=0.641, a(45)=0.676 (journal, L=0.3) ==")
    print(df.round(3).to_string(index=False)); print(params)
    OUTJ["A3"] = dict(table=df.to_dict(orient="records"), params=params)
    df.to_csv(OUT / "A3_acceptance_forms.csv", index=False)
    # curvature lambda_a with the level anchored at the 28-day arm
    rows = []
    for lam in (0.05, 0.075, LAM0, 0.15, 0.20):
        a = acc_exp(lam)
        for delta in (0.0025, 0.005, 0.01):
            V = value(delta, 0.3, DA0, DL0, 0.0, a)[0]
            rows.append(dict(lam_a=lam, abar=float(0.641 / (1 - math.exp(-28 * lam))), delta=delta, T_star=int(TS[np.argmax(V)]), a14=float(a[0][TS == 14][0])))
    d2 = pd.DataFrame(rows)
    print(d2.round(3).to_string(index=False))
    OUTJ["A3_lambda"] = d2.to_dict(orient="records")


# ========================================================= A4 pass-through ==
def closed_form(q=1.0, L=0.3, delta=0.005, lam=0.099, m=1.0):
    Ad = 1 + lam * (q + L + delta * m) / delta
    return float((Ad - lambertw(np.exp(Ad)).real) / lam)


def a4():
    # (i) pass-through implied by the RCT medians: 48 d at T=45, 36 d at T=28
    rho_rct = (48 - 36) / (45 - 28)
    # (ii) generalised closed form: E tau = tau0 + rho (T - T0)  <=>  delta' = rho*delta, m' = T0 - tau0/rho
    rows = []
    for rho in (1.0, 0.93, rho_rct):
        for delta in (0.0025, 0.005, 0.01):
            for L in (0.1, 0.3, 0.6):
                tau0, T0 = 36.0, 28.0
                dp, mp = rho * delta, T0 - tau0 / rho
                Tc = closed_form(L=L, delta=dp, lam=LAM0, m=mp)
                Tg = np.linspace(1, 200, 400000)
                a = ABAR0 * (1 - np.exp(-LAM0 * Tg))
                Vn = a * (1 - delta * (tau0 + rho * (Tg - T0))) - (1 - a) * L
                Tn = float(Tg[np.argmax(Vn)])
                rows.append(dict(rho=rho, delta=delta, L=L, T_closed=Tc, T_numeric=Tn, diff=abs(Tc - Tn)))
    df = pd.DataFrame(rows)
    print("\n== A4 partial anchoring: journal delivery anchored on the RCT (median 36 d at T=28), Delta=0 ==")
    print(df.pivot_table(index=["rho", "L"], columns="delta", values="T_closed").round(1).to_string())
    print("max |closed - numeric| =", df["diff"].max(), "| rho_rct =", rho_rct)
    # (iii) the same pass-through applied to the ICLR profile (with the length penalty)
    r2 = []
    a = (ABAR0 * (1 - np.exp(-LAM0 * TS)))[None, :]
    for rho in (1.0, 0.93, 0.85, rho_rct):
        for delta in (0.0025, 0.005, 0.01):
            for eta in (0.0, 1.3):
                V = value(delta, 0.3, DA0, DL0, eta, a, rho=rho)[0]
                r2.append(dict(rho=rho, delta=delta, eta=eta, T_star=int(TS[np.argmax(V)])))
    d2 = pd.DataFrame(r2)
    print(d2.pivot_table(index=["rho", "eta"], columns="delta", values="T_star").to_string())
    OUTJ["A4"] = dict(rho_rct=rho_rct, closed=df.to_dict(orient="records"), max_diff=float(df["diff"].max()), profile=d2.to_dict(orient="records"))
    df.to_csv(OUT / "A4_passthrough_closed.csv", index=False); d2.to_csv(OUT / "A4_passthrough_profile.csv", index=False)


# ================================================================ A5 tails ==
def a5():
    rows = []
    aJ = (ABAR0 * (1 - np.exp(-LAM0 * TS)))[None, :]
    one = np.ones((1, len(TS)))
    for tail in TAILS:
        atom, late, spre, reach = PRE[tail]
        big = pre_moments(tail, ts=np.array([400]))
        r = dict(tail=tail, reach48_T12=float(reach[TS == 12][0]), reach48_T21=float(reach[TS == 21][0]), reach48_T45=float(reach[TS == 45][0]),
                 reach48_T90=float(reach[TS == 90][0]), reach48_T400=float(big[3][0]),
                 atom_T21=float(atom[TS == 21][0]), late_T21=float(late[TS == 21][0]), atom_T45=float(atom[TS == 45][0]), late_T45=float(late[TS == 45][0]))
        for name, acc, L, eta in [("conf_eta1.3", one, 0.0, 1.3), ("conf_eta0", one, 0.0, 0.0), ("jour_eta0", aJ, 0.3, 0.0), ("jour_eta1.3", aJ, 0.3, 1.3)]:
            for delta in (0.0025, 0.005, 0.01):
                V = value(delta, L, DA0, DL0, eta, acc, tail=tail)[0]
                r[f"T_{name}_{delta}"] = int(TS[np.argmax(V)])
        rows.append(r)
    df = pd.DataFrame(rows)
    print("\n== A5 tail of the hazard profile ==")
    print(df.round(3).T.to_string()); print("exp tail: h = exp(%.3f - %.3f s)" % (EXP_C, EXP_K))
    OUTJ["A5"] = dict(table=df.to_dict(orient="records"), exp_c=float(EXP_C), exp_k=float(EXP_K), pl_a=PL_A, pl_b=PL_B)
    df.to_csv(OUT / "A5_tails.csv", index=False)


# ================================================================ A6 Delta ==
def a6():
    rows = []
    aJ = (ABAR0 * (1 - np.exp(-LAM0 * TS)))[None, :]
    one = np.ones((1, len(TS)))
    for D in (0.0, 0.02, 0.05, 0.085, 0.15, 0.30):
        for eta in (0.0, 1.3):
            for delta in (0.0025, 0.005, 0.01):
                Vj = value(delta, 0.3, D, D, eta, aJ)[0]; Vc = value(delta, 0.0, D, D, eta, one)[0]
                rows.append(dict(Delta=D, eta=eta, delta=delta, T_journal=int(TS[np.argmax(Vj)]), T_conference=int(TS[np.argmax(Vc)])))
    df = pd.DataFrame(rows)
    print("\n== A6 compression penalty from 0 ==")
    print(df.pivot_table(index=["Delta", "eta"], columns="delta", values=["T_journal", "T_conference"]).to_string())
    OUTJ["A6"] = df.to_dict(orient="records")
    df.to_csv(OUT / "A6_delta.csv", index=False)


# ===================================================== A7/A8 joint analysis ==
def lhs(n, k, g):
    u = (np.argsort(g.random((n, k)), axis=0) + g.random((n, k))) / n
    return u


def first_order(x, y, bins=20):
    """First-order sensitivity index Var(E[y|x]) / Var(y), x continuous (equal-count bins) or categorical."""
    if x.dtype.kind in "iu" and len(np.unique(x)) <= 10:
        groups = [y[x == v] for v in np.unique(x)]
    else:
        q = np.quantile(x, np.linspace(0, 1, bins + 1)); q[-1] += 1e-12
        idx = np.clip(np.searchsorted(q, x, side="right") - 1, 0, bins - 1)
        groups = [y[idx == b] for b in range(bins)]
    m = np.array([g.mean() for g in groups]); w = np.array([len(g) for g in groups]) / len(y)
    return float((w * (m - y.mean()) ** 2).sum() / y.var())


def a7(n=20000):
    g = np.random.Generator(np.random.PCG64(SEED))
    u = lhs(n, 8, g)
    delta = 0.0025 * (0.01 / 0.0025) ** u[:, 0]            # log-uniform 0.0025-0.01
    L = 0.1 + 0.5 * u[:, 1]                                 # 0.1-0.6
    Delta = 0.30 * u[:, 2]                                  # 0-0.30
    eta = 2.0 * u[:, 3]                                     # 0-2
    lam = 0.05 * (0.20 / 0.05) ** u[:, 4]                   # log-uniform 0.05-0.20, level anchored at a(28)
    rho = 0.7 + 0.3 * u[:, 5]                               # pass-through 0.7-1
    mu = 0.20 + 0.22 * u[:, 6]                              # late hazard 0.20-0.42 (range of the cohorts)
    tail = np.minimum((u[:, 7] * 4).astype(int), 3)
    res = {}
    Tc = np.zeros(n, int); Tj = np.zeros(n, int); uni_c = np.zeros(n, bool); uni_j = np.zeros(n, bool)
    lossJ14 = np.zeros(n); lossJ42 = np.zeros(n); lossC14 = np.zeros(n); lossC42 = np.zeros(n)
    i14, i42 = int(np.where(TS == 14)[0][0]), int(np.where(TS == 42)[0][0])
    for ti, tname in enumerate(TAILS):
        sel = tail == ti
        one = np.ones((sel.sum(), len(TS)))
        Vc = value(delta[sel], 0.0, Delta[sel], Delta[sel], eta[sel], one, tail=tname, mu=mu[sel], rho=1.0)
        Vj = value(delta[sel], L[sel], Delta[sel], Delta[sel], eta[sel], acc_exp(lam[sel]), tail=tname, mu=mu[sel], rho=rho[sel])
        Tc[sel] = tstar(Vc); Tj[sel] = tstar(Vj)
        for V, uni in ((Vc, uni_c), (Vj, uni_j)):
            d = np.sign(np.round(np.diff(V, axis=1), 12))
            # unimodal: once the increments turn negative they never turn positive again
            first_neg = np.argmax(d < 0, axis=1)
            has_neg = (d < 0).any(axis=1)
            later_pos = np.array([(row[fn:] > 0).any() if hn else False for row, fn, hn in zip(d, first_neg, has_neg)])
            uni[sel] = ~later_pos
        lossJ14[sel] = 100 * (Vj.max(1) - Vj[:, i14]); lossJ42[sel] = 100 * (Vj.max(1) - Vj[:, i42])
        lossC14[sel] = 100 * (Vc.max(1) - Vc[:, i14]); lossC42[sel] = 100 * (Vc.max(1) - Vc[:, i42])
    qs = [5, 25, 50, 75, 95]
    summ = dict(n=n,
                conference=dict(q={q: float(np.percentile(Tc, q)) for q in qs}, share_corner3=float((Tc == 3).mean()), share_le7=float((Tc <= 7).mean()),
                                share_14_25=float(((Tc >= 14) & (Tc <= 25)).mean()), share_at90=float((Tc == 90).mean()), unimodal=float(uni_c.mean()),
                                loss14_gt_loss42=float((lossC14 > lossC42).mean())),
                journal=dict(q={q: float(np.percentile(Tj, q)) for q in qs}, share_25_44=float(((Tj >= 25) & (Tj <= 44)).mean()), share_lt21=float((Tj < 21).mean()),
                             share_gt60=float((Tj > 60).mean()), share_at90=float((Tj == 90).mean()), unimodal=float(uni_j.mean()),
                             loss14_gt_loss42=float((lossJ14 > lossJ42).mean()),
                             loss14_q={q: float(np.percentile(lossJ14, q)) for q in qs}, loss42_q={q: float(np.percentile(lossJ42, q)) for q in qs}))
    # conference with and without compression growth
    for lab, m in (("eta<0.5", eta < 0.5), ("eta>=1", eta >= 1.0)):
        summ["conference"][lab] = {q: float(np.percentile(Tc[m], q)) for q in qs}
    Xc = dict(delta=delta, Delta=Delta, eta=eta, mu=mu, tail=tail)
    Xj = dict(delta=delta, L=L, Delta=Delta, eta=eta, lam_a=lam, rho=rho, mu=mu, tail=tail)
    S_c = {k: first_order(v, Tc.astype(float)) for k, v in Xc.items()}
    S_j = {k: first_order(v, Tj.astype(float)) for k, v in Xj.items()}
    S_c["interactions"] = 1 - sum(S_c.values()); S_j["interactions"] = 1 - sum(S_j.values())
    summ["S_conference"] = S_c; summ["S_journal"] = S_j
    # conditional medians of the journal optimum by lambda_a tercile and delta tercile (map)
    lt = np.quantile(lam, [1 / 3, 2 / 3]); dt = np.quantile(delta, [1 / 3, 2 / 3])
    mp = {}
    for i, ln in enumerate(["low", "mid", "high"]):
        for j, dn in enumerate(["low", "mid", "high"]):
            m = (np.digitize(lam, lt) == i) & (np.digitize(delta, dt) == j)
            mp[f"lam_{ln}|delta_{dn}"] = [float(np.percentile(Tj[m], q)) for q in (5, 50, 95)]
    summ["journal_map"] = mp; summ["lam_terciles"] = [float(x) for x in lt]; summ["delta_terciles"] = [float(x) for x in dt]
    print("\n== A7/A8 joint sensitivity (LHS, n=%d) ==" % n)
    print(json.dumps(summ, indent=1))
    OUTJ["A7"] = summ
    np.savez_compressed(OUT / "A7_draws.npz", delta=delta, L=L, Delta=Delta, eta=eta, lam=lam, rho=rho, mu=mu, tail=tail, Tc=Tc, Tj=Tj)
    return dict(Tc=Tc, Tj=Tj, S_c=S_c, S_j=S_j, eta=eta)


# ============================================================ A9 E4 extras ==
def a9():
    by = pd.DataFrame(RES["E4"]["by_year"])
    out = {}
    for outcome in ("log_words", "confidence", "score", "dev_panel"):
        d = by[(by["outcome"] == outcome) & (by["bin"] == "last24h")].dropna(subset=["coef", "se"])
        def pool(years):
            x = d[d["year"].isin(years)]
            w = 1 / x["se"] ** 2
            est = float((w * x["coef"]).sum() / w.sum()); se = float(math.sqrt(1 / w.sum()))
            Q = float((w * (x["coef"] - est) ** 2).sum()); k = len(x)
            return dict(est=est, se=se, k=k, Q=Q, p_het=float(stats.chi2.sf(Q, k - 1)) if k > 1 else float("nan"))
        pre, post = pool(range(2018, 2024)), pool([2024, 2025, 2026])
        post_no25 = pool([2024, 2026])
        z = (post["est"] - pre["est"]) / math.sqrt(pre["se"] ** 2 + post["se"] ** 2)
        out[outcome] = dict(pre=pre, post=post, post_excl_2025=post_no25, diff=post["est"] - pre["est"], z=float(z), p=float(2 * stats.norm.sf(abs(z))),
                            by_year={int(r.year): [float(r.coef), float(r.se)] for r in d.itertuples()})
    # same-form contrasts for length: 2022 vs 2023 (pre-ChatGPT, forms differ), 2024-2026 share one form
    pooled = pd.DataFrame(RES["E4"]["pooled"])
    r = pooled[(pooled["outcome"] == "dev_panel") & (pooled["bin"] == "last24h")].iloc[0]
    coef, se, mean_early = float(r["coef"]), float(r["se"]), float(r["mean_early"])
    ci90 = (coef - 1.6449 * se, coef + 1.6449 * se); ci95 = (coef - 1.96 * se, coef + 1.96 * se)
    tost = {}
    for margin in (0.03, 0.05, 0.10):
        p_low = stats.norm.sf((coef + margin) / se); p_high = stats.norm.cdf((coef - margin) / se)
        tost[str(margin)] = dict(p=float(max(p_low, p_high)), margin_pct_of_mean=100 * margin / mean_early)
    out["dev_panel_equivalence"] = dict(coef=coef, se=se, mean_early=mean_early, ci90=ci90, ci95=ci95, tost=tost)
    for b in ("mid", "last48h", "late"):
        rr = pooled[(pooled["outcome"] == "dev_panel") & (pooled["bin"] == b)].iloc[0]
        out["dev_panel_equivalence"][b] = [float(rr["coef"]), float(rr["se"])]
    print("\n== A9 timing penalties before/after ChatGPT; equivalence of panel deviation ==")
    print(json.dumps(out, indent=1))
    OUTJ["A9"] = out


# ============================================================== A10 kappa ===
def a10():
    k = pd.DataFrame(RES["E6"]["kappa"])
    acc = k[k["pair"] == "all"].set_index("year")["accept_rate"]
    rows = []
    for y in range(2018, 2026):
        d = k[k["year"] == y].set_index("pair")
        p = float(acc[y]); pe = p * p + (1 - p) ** 2
        def se(kap, n):
            po = pe + kap * (1 - pe)
            return math.sqrt(po * (1 - po) / (n * (1 - pe) ** 2))
        r = dict(year=y)
        for pair in ("early-early", "atom-atom", "late-late", "early-late", "early-atom", "atom-late", "all"):
            r[pair] = float(d.loc[pair, "kappa"]); r["n_" + pair] = int(d.loc[pair, "n_pairs"]); r["se_" + pair] = se(r[pair], r["n_" + pair])
        rows.append(r)
    df = pd.DataFrame(rows)
    res = {}
    for yrs, lab in ((range(2021, 2026), "2021-2025"), (range(2018, 2026), "2018-2025")):
        x = df[df["year"].isin(yrs)]
        for a, b in (("early-early", "late-late"), ("early-early", "atom-atom"), ("atom-atom", "late-late")):
            dif = x[a] - x[b]; v = x["se_" + a] ** 2 + x["se_" + b] ** 2
            w = 1 / v
            est = float((w * dif).sum() / w.sum()); se_ = float(math.sqrt(1 / w.sum()))
            res[f"{a} minus {b} | {lab}"] = dict(est=est, se=se_, z=est / se_, p=float(2 * stats.norm.sf(abs(est / se_))),
                                                 n_positive=int((dif > 0).sum()), k=len(x), by_year={int(yy): float(dd) for yy, dd in zip(x["year"], dif)})
    print("\n== A10 agreement by timing pair (approximate SE: pairs treated as independent) ==")
    print(df[["year", "early-early", "n_early-early", "se_early-early", "atom-atom", "se_atom-atom", "late-late", "n_late-late", "se_late-late"]].round(3).to_string(index=False))
    print(json.dumps(res, indent=1))
    OUTJ["A10"] = dict(table=df.to_dict(orient="records"), pooled=res)


# ============================================================ A11 completion =
def a11():
    with gzip.open(REPO / "data" / "processed" / "papers.csv.gz") as f:
        pap = pd.read_csv(io.TextIOWrapper(f, encoding="utf8"))
    pap = pap[pap["n_reviews"] >= 1]
    rows = []
    for y, g in pap.groupby("year"):
        n = g["n_reviews"]
        rows.append(dict(year=int(y), manuscripts=len(g), mean_reviews=float(n.mean()), modal=int(n.mode().iloc[0]), share_lt3=float((n < 3).mean()),
                         share_eq1=float((n == 1).mean()), share_ge4=float((n >= 4).mean())))
    df = pd.DataFrame(rows)
    print("\n== A11 reviews per manuscript ==")
    print(df.round(3).to_string(index=False))
    OUTJ["A11"] = df.to_dict(orient="records")
    df.to_csv(OUT / "A11_reviews_per_manuscript.csv", index=False)


def figure(d):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9, "legend.fontsize": 8, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5})
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.5))
    ax = axes[0]
    b = np.arange(2.5, 91.5, 1)
    ax.hist(d["Tc"][d["eta"] < 0.5], bins=b, alpha=0.75, label=r"conference, $\eta<0.5$", density=True)
    ax.hist(d["Tc"][d["eta"] >= 1.0], bins=b, alpha=0.75, label=r"conference, $\eta\geq1$", density=True)
    ax.set_xlim(0, 62); ax.set_xlabel("Optimal window T* (days)"); ax.set_ylabel("Density over joint draws")
    ax.set_title("(a) Conference: T* across joint draws"); ax.legend(frameon=False)
    ax = axes[1]
    ax.hist(d["Tj"], bins=b, alpha=0.8, color="C2", density=True, label="journal")
    for q, ls in ((5, ":"), (50, "-"), (95, ":")):
        ax.axvline(np.percentile(d["Tj"], q), color="k", lw=0.9, ls=ls)
    ax.axvspan(23, 44, color="grey", alpha=0.15, lw=0, label="scenario table,\nPanel B (23\u201344)")
    ax.set_xlim(0, 91); ax.set_xlabel("Optimal window T* (days)"); ax.set_ylabel("Density over joint draws")
    ax.set_title("(b) Journal, exponential acceptance curves"); ax.legend(frameon=False)
    ax = axes[2]
    names = {"delta": r"delay cost $\delta$", "L": "decline loss $L$", "Delta": r"shortfall $\Delta$", "eta": r"elasticity $\eta$", "lam_a": r"acceptance curvature $\lambda_a$",
             "rho": r"pass-through $\rho$", "mu": r"late hazard $\mu$", "tail": "tail of the hazard", "interactions": "interactions"}
    keys = ["lam_a", "delta", "eta", "Delta", "L", "rho", "tail", "mu", "interactions"]
    y = np.arange(len(keys))
    ax.barh(y - 0.2, [d["S_j"].get(k, 0) for k in keys], height=0.38, color="C2", label="journal")
    ax.barh(y + 0.2, [d["S_c"].get(k, 0) for k in keys], height=0.38, color="C0", label="conference")
    ax.set_yticks(y); ax.set_yticklabels([names[k] for k in keys], fontsize=8); ax.invert_yaxis()
    ax.set_xlabel("Share of the variance of T*"); ax.set_title("(c) What moves the optimum"); ax.legend(frameon=False)
    for a_ in axes:
        for sp in ("top", "right"): a_.spines[sp].set_visible(False)
    fig.tight_layout(); fig.savefig(OUT / "fig6_joint_sensitivity.pdf"); fig.savefig(OUT / "fig6_joint_sensitivity.png", dpi=200); plt.close(fig)




# ================================================================ A12 extras =
def a12():
    out = {}
    out["closed_form_lam0"] = closed_form(L=0.3, delta=0.005, lam=LAM0, m=1.0)
    out["closed_form_lam010"] = closed_form(L=0.3, delta=0.005, lam=0.10, m=1.0)
    by = pd.DataFrame(RES["E4"]["by_year"])
    d = by[(by["outcome"] == "log_words") & (by["bin"] == "last24h")]
    def pool(years):
        x = d[d["year"].isin(years)]; w = 1 / x["se"] ** 2
        return float((w * x["coef"]).sum() / w.sum()), float(math.sqrt(1 / w.sum()))
    pre_no23 = pool(range(2018, 2023)); post = pool([2024, 2025, 2026])
    out["log_words_pre_excl_2023"] = dict(est=pre_no23[0], se=pre_no23[1], z_vs_post=(post[0] - pre_no23[0]) / math.sqrt(pre_no23[1] ** 2 + post[1] ** 2))
    a1t = {r["stat"]: r for r in OUTJ["A1"]["table"]}
    pt = -(a1t["share_atom"]["slope"] + a1t["share_late"]["slope"]) * DA0
    ub = -(a1t["share_atom"]["lo"] + a1t["share_late"]["lo"]) * DA0
    out["marginal_quality_per_day_pct"] = dict(point=100 * pt, upper=100 * ub)
    out["sign_test_5of5_two_sided"] = 2 * 0.5 ** 5
    out["passthrough_median_lower_bound"] = 1 - a1t["s_q50"]["hi"]
    # share of the window's extra day that shows up as delay at the median: 1 - slope
    out["passthrough_median"] = [1 - a1t["s_q50"]["hi"], 1 - a1t["s_q50"]["slope"], 1 - a1t["s_q50"]["lo"]]
    # ranges of T* in the submitted table by eta
    ts_ = pd.DataFrame(RES["E7"]["tstar"])
    for conf in (True, False):
        for eta in (0.0, 1.3):
            x = ts_[(ts_["conference"] == conf) & (ts_["eta"] == eta)]
            out[f"table_range_conf{conf}_eta{eta}"] = [int(x["T_star"].min()), int(x["T_star"].max())]
            for T in (7, 14, 21, 28, 42):
                out[f"loss{T}_conf{conf}_eta{eta}"] = [float(x[f"loss_T{T}"].min()), float(x[f"loss_T{T}"].max())]
    print("\n== A12 extras =="); print(json.dumps(out, indent=1))
    OUTJ["A12"] = out


if __name__ == "__main__":
    a1(); a2(); a3(); a4(); a5(); a6()
    d = a7()
    a9(); a10(); a11(); a12()
    figure(d)
    (OUT / "revision_results.json").write_text(json.dumps(OUTJ, indent=1, default=float))
    print("\nwritten", OUT / "revision_results.json")
