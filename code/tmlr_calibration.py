"""
Journal calibration with a delivery profile measured in a journal (TMLR) instead of the conference profile.

usage: python tmlr_calibration.py <repository root> <tmlr output folder>

Reads output/results.json (conference profile, penalties, submitted scenarios) and tmlr_results.json (hazard by day
relative to the deadline in the same-day panels of TMLR). Writes tmlr_calibration.json and a CSV.
The acceptance curve is the one of the manuscript (two arms of Chetty et al. 2014).
"""
import json, sys, math, csv, os
import numpy as np
from scipy.optimize import brentq

REPO, OUT = sys.argv[1], sys.argv[2]
RES = json.load(open(os.path.join(REPO, "output", "results.json")))
TM = json.load(open(os.path.join(OUT, "tmlr_results.json")))
TS = np.arange(3, 91); MAX_LATE, T_REF = 60, 20.5
PROF = RES["E7"]["profile"]; H_ICLR = {int(k): float(v) for k, v in PROF["h"].items()}
MU_ICLR, PL_A, PL_B = float(PROF["mu"]), float(PROF["pl_a"]), float(PROF["pl_b"])
DA0, DL0 = float(RES["E7"]["Delta_atom"]), float(RES["E7"]["Delta_late"])

def h_iclr(s):
    return H_ICLR[s] if s < 15 else min(0.999, PL_A * (1 + s) ** (-PL_B))

# TMLR: pooled 14- and 28-day panels for the last 14 days, 28-day panels for days 14-27, constant beyond the data
hp = {x["day"]: x for x in TM["hazard"]["pooled_last14"] if x["day"] >= 0}
h28 = {x["day"]: x for x in TM["hazard"]["W28"] if 14 <= x["day"] <= 27}
FLOOR = sum(x["events"] for x in h28.values()) / sum(x["at_risk"] for x in h28.values())
MU_TMLR = TM["late_hazard"]["pooled"]["1-60"]["hazard"]
def h_tmlr(s):
    if s <= 13: return hp[s]["hazard"]
    if s <= 27: return h28[s]["hazard"]
    return FLOOR
def h_tmlr_zero_beyond(s):
    return h_tmlr(s) if s <= 27 else 0.0

def pre(hfun):
    atom, late, spre = [], [], []
    for T in TS:
        surv = 1.0; a = sp = 0.0
        for s in range(int(T) - 1, -1, -1):
            h = hfun(s); p = surv * h
            if s < 2: a += p
            sp += s * p; surv *= 1 - h
        atom.append(a); late.append(surv); spre.append(sp)
    return np.array(atom), np.array(late), np.array(spre)

def mean_days_late(mu):
    d = np.arange(1, MAX_LATE)
    return float((d * mu * (1 - mu) ** (d - 1)).sum() + MAX_LATE * (1 - mu) ** (MAX_LATE - 1))

def profile(hfun, mu):
    atom, late, spre = pre(hfun)
    delay = TS - (spre - late * mean_days_late(mu) + 0.5)
    return dict(atom=atom, late=late, delay=delay)

def value(P, delta, L, Da, Dl, eta, acc):
    comp = (T_REF / TS) ** eta
    Eq = 1.0 - Da * comp * P["atom"] - Dl * comp * P["late"]
    return acc * (Eq - delta * P["delay"]) - (1 - acc) * L

a28, a45 = 0.641, 0.676
lam = brentq(lambda l: (1 - math.exp(-45 * l)) / (1 - math.exp(-28 * l)) - a45 / a28, 1e-4, 1.0)
abar = a28 / (1 - math.exp(-28 * lam)); ACC = abar * (1 - np.exp(-lam * TS))

P_I = profile(h_iclr, MU_ICLR); P_T = profile(h_tmlr, MU_TMLR); P_T0 = profile(h_tmlr_zero_beyond, MU_TMLR)
DL_T = 1 - math.exp(TM["fe_log_chars_late_vs_on_time"]["b"])       # late reports shorter than on-time ones (same manuscript)
sel = TS <= 60
rows = []
for eta in (0.0, 1.3):
    for L in (0.1, 0.3, 0.6):
        for delta in (0.0025, 0.005, 0.01):
            t = lambda V: int(TS[sel][np.argmax(V[sel])])
            r = dict(eta=eta, L=L, delta=delta,
                     conference_profile=t(value(P_I, delta, L, DA0, DL0, eta, ACC)),
                     journal_timing_conference_penalties=t(value(P_T, delta, L, DA0, DL0, eta, ACC)),
                     journal_timing_journal_penalties=t(value(P_T, delta, L, 0.0, DL_T, eta, ACC)),
                     journal_timing_no_penalty=t(value(P_T, delta, L, 0.0, 0.0, eta, ACC)),
                     journal_timing_zero_hazard_beyond_27=t(value(P_T0, delta, L, 0.0, DL_T, eta, ACC)))
            rows.append(r)
i = lambda T: int(np.where(TS == T)[0][0])
out = dict(lam_a=lam, abar=abar, mu_tmlr=MU_TMLR, mu_iclr=MU_ICLR, hazard_floor_beyond_27=FLOOR, Delta_late_tmlr=DL_T,
           expected_delay={k: {str(T): float(P["delay"][i(T)]) for T in (14, 21, 28, 35, 42)} for k, P in (("conference", P_I), ("journal", P_T))},
           late_share={k: {str(T): float(P["late"][i(T)]) for T in (14, 21, 28, 35, 42)} for k, P in (("conference", P_I), ("journal", P_T))},
           final48h_share={k: {str(T): float(P["atom"][i(T)]) for T in (14, 21, 28, 35, 42)} for k, P in (("conference", P_I), ("journal", P_T))},
           slope_delay_28_42={k: float((P["delay"][i(42)] - P["delay"][i(28)]) / 14) for k, P in (("conference", P_I), ("journal", P_T))},
           scenarios=rows)
d = np.array([[r["journal_timing_journal_penalties"] - r["conference_profile"] for r in rows]])
out["max_abs_shift_days"] = int(np.abs(d).max()); out["shift_range"] = [int(d.min()), int(d.max())]
d0 = np.array([r["journal_timing_zero_hazard_beyond_27"] - r["conference_profile"] for r in rows])
out["shift_range_zero_hazard_beyond_27"] = [int(d0.min()), int(d0.max())]
out["range_journal_profile_zero_hazard_beyond_27"] = [min(r["journal_timing_zero_hazard_beyond_27"] for r in rows), max(r["journal_timing_zero_hazard_beyond_27"] for r in rows)]
out["range_conference_profile"] = [min(r["conference_profile"] for r in rows), max(r["conference_profile"] for r in rows)]
out["range_journal_profile"] = [min(r["journal_timing_journal_penalties"] for r in rows), max(r["journal_timing_journal_penalties"] for r in rows)]
# model check: the profile built from the pooled hazards reproduces the observed shares of the two windows
out["check_profile_vs_observed"] = {"W14": dict(model_late=float(P_T["late"][i(14)]), observed_late=TM["groups"]["W14"]["share_late"], model_final48h=float(P_T["atom"][i(14)]), observed_final48h=TM["groups"]["W14"]["share_final48h"], model_delay=float(P_T["delay"][i(14)]), observed_mean_tau=TM["groups"]["W14"]["tau_mean"]),
                                    "W28": dict(model_late=float(P_T["late"][i(28)]), observed_late=TM["groups"]["W28"]["share_late"], model_final48h=float(P_T["atom"][i(28)]), observed_final48h=TM["groups"]["W28"]["share_final48h"], model_delay=float(P_T["delay"][i(28)]), observed_mean_tau=TM["groups"]["W28"]["tau_mean"])}
json.dump(out, open(os.path.join(OUT, "tmlr_calibration.json"), "w"), indent=1)
with open(os.path.join(OUT, "tmlr_calibration_scenarios.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
for r in rows: print(r)
print({k: v for k, v in out.items() if k != "scenarios"})
