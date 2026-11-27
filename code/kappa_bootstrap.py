"""Agreement by timing pair with a manuscript-clustered bootstrap (replaces the approximate standard errors of A10).

usage: python code/kappa_bootstrap.py <repository root> <reviews.csv or reviews.csv.gz> <output folder>

Cohen's kappa on the accept vote (score >= 6) between reviews of the same manuscript, by the timing of the two
reviews (early: four days or more before the deadline; final 48 h; late), as in exp_quality.e6. Every kappa is a
function of sums over manuscripts, so a bootstrap resample of manuscripts is a reweighting of those sums.
The script first reproduces the archived kappas of output/results.json. Seed PCG64(20261127); 4,000 resamples.
"""
import json, sys, os, math
import numpy as np
import pandas as pd

REPO, REV, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(OUT, exist_ok=True)
RES = json.load(open(os.path.join(REPO, "output", "results.json")))
SEED, B = 20261127, 4000
GROUPS = ["early", "atom", "late"]
PAIRS = [("early", "early"), ("atom", "atom"), ("late", "late"), ("early", "atom"), ("early", "late"), ("atom", "late")]

df = pd.read_csv(REV, low_memory=False, usecols=["year", "forum", "score", "bin", "in_sample", "n_reviews"])
df = df[df["in_sample"].astype(bool) & (df["n_reviews"] >= 2) & df["score"].notna()].copy()
df["accept"] = (df["score"] >= 6.0).astype(int)
df["grp"] = df["bin"].map({"early": "early", "mid": "mid", "last48h": "atom", "last24h": "atom", "late": "late"})

def components(g):
    """per-manuscript sums from which every kappa is computed; columns indexed by name"""
    forums = pd.Index(g["forum"].unique())
    n = {}; K = {}
    for k in GROUPS + ["all"]:
        gg = g if k == "all" else g[g["grp"] == k]
        a = gg.groupby("forum")["accept"].agg(["sum", "size"]).reindex(forums).fillna(0.0)
        K[k] = a["sum"].to_numpy(float); n[k] = a["size"].to_numpy(float)
    X = {}
    for k in GROUPS + ["all"]:
        ok = n[k] >= 2
        X[f"{k}|agree"] = np.where(ok, K[k] * (K[k] - 1) / 2 + (n[k] - K[k]) * (n[k] - K[k] - 1) / 2, 0.0)
        X[f"{k}|pairs"] = np.where(ok, n[k] * (n[k] - 1) / 2, 0.0)
        X[f"{k}|c"] = np.where(ok, K[k] * (n[k] - 1), 0.0)
        X[f"{k}|d"] = np.where(ok, n[k] * (n[k] - 1), 0.0)
    for a, b in PAIRS:
        if a == b:
            continue
        X[f"{a}-{b}|agree"] = K[a] * K[b] + (n[a] - K[a]) * (n[b] - K[b])
        X[f"{a}-{b}|pairs"] = n[a] * n[b]
        X[f"{a}-{b}|p1"] = K[a] * n[b]
        X[f"{a}-{b}|p2"] = K[b] * n[a]
    names = list(X); M = np.column_stack([X[c] for c in names])
    return names, M

def kappas(S, names):
    """S: array (..., ncols) of (weighted) sums -> dict pair -> kappa array, and n_pairs"""
    ix = {c: i for i, c in enumerate(names)}; out = {}; npairs = {}
    for a, b in PAIRS + [("all", "all")]:
        if a == b:
            po = S[..., ix[f"{a}|agree"]] / S[..., ix[f"{a}|pairs"]]; p1 = p2 = S[..., ix[f"{a}|c"]] / S[..., ix[f"{a}|d"]]
            key = "all" if a == "all" else f"{a}-{b}"; npairs[key] = S[..., ix[f"{a}|pairs"]]
        else:
            key = f"{a}-{b}"; pr = S[..., ix[f"{key}|pairs"]]
            po = S[..., ix[f"{key}|agree"]] / pr; p1 = S[..., ix[f"{key}|p1"]] / pr; p2 = S[..., ix[f"{key}|p2"]] / pr; npairs[key] = pr
        pe = p1 * p2 + (1 - p1) * (1 - p2)
        out[key] = (po - pe) / (1 - pe)
    return out, npairs

rng = np.random.default_rng(SEED)
arch = {(r["year"], r["pair"]): r for r in RES["E6"]["kappa"]}
est, reps, rows, maxdiff = {}, {}, [], 0.0
for y in range(2018, 2026):
    g = df[df["year"] == y]
    names, M = components(g)
    k0, n0 = kappas(M.sum(axis=0), names)
    for key in k0:
        a = arch[(y, key)]
        maxdiff = max(maxdiff, abs(float(k0[key]) - a["kappa"])); assert int(round(float(n0[key]))) == a["n_pairs"], (y, key)
    W = rng.multinomial(M.shape[0], np.full(M.shape[0], 1.0 / M.shape[0]), size=B).astype(float)
    kb, _ = kappas(W @ M, names)
    est[y] = {k: float(v) for k, v in k0.items()}; reps[y] = kb
    r = dict(year=y, manuscripts=int(M.shape[0]))
    for key in ("early-early", "atom-atom", "late-late", "all"):
        r[key] = est[y][key]; r["n_" + key] = int(round(float(n0[key]))); r["se_" + key] = float(np.nanstd(kb[key], ddof=1))
    for a, b in (("early-early", "late-late"), ("early-early", "atom-atom"), ("atom-atom", "late-late")):
        d = kb[a] - kb[b]; nm = f"{a} minus {b}"
        r[nm] = est[y][a] - est[y][b]; r["se_" + nm] = float(np.nanstd(d, ddof=1)); r["lo_" + nm], r["hi_" + nm] = [float(x) for x in np.nanpercentile(d, [2.5, 97.5])]
    rows.append(r)
tab = pd.DataFrame(rows)
pooled = {}
for yrs, lab in ((range(2021, 2026), "2021-2025"), (range(2018, 2026), "2018-2025")):
    x = tab[tab["year"].isin(yrs)]
    for a, b in (("early-early", "late-late"), ("early-early", "atom-atom"), ("atom-atom", "late-late")):
        nm = f"{a} minus {b}"
        w = 1 / x["se_" + nm].to_numpy() ** 2; w = w / w.sum()
        e = float((w * x[nm].to_numpy()).sum())
        rb = sum(wi * (reps[int(y)][a] - reps[int(y)][b]) for wi, y in zip(w, x["year"]))      # cohorts resampled independently
        se = float(np.nanstd(rb, ddof=1)); lo, hi = [float(v) for v in np.nanpercentile(rb, [2.5, 97.5])]
        # unweighted mean of the cohort differences, as a check that the weights do not drive the result
        ru = sum((reps[int(y)][a] - reps[int(y)][b]) / len(x) for y in x["year"])
        pooled[f"{nm} | {lab}"] = dict(est=e, se=se, lo=lo, hi=hi, z=e / se, share_boot_positive=float(np.mean(rb > 0)),
                                       n_positive=int((x[nm] > 0).sum()), k=int(len(x)),
                                       unweighted=float(x[nm].mean()), unweighted_lo=float(np.nanpercentile(ru, 2.5)), unweighted_hi=float(np.nanpercentile(ru, 97.5)),
                                       by_year={int(y): float(v) for y, v in zip(x["year"], x[nm])})
out = dict(seed=SEED, n_boot=B, max_abs_diff_with_archived_kappa=maxdiff, table=tab.to_dict(orient="records"), pooled=pooled)
json.dump(out, open(os.path.join(OUT, "kappa_bootstrap.json"), "w"), indent=1)
tab.to_csv(os.path.join(OUT, "A13_kappa_bootstrap.csv"), index=False)
print("reproduction of the archived kappas: max abs difference", maxdiff)
print(tab[["year", "manuscripts", "early-early", "se_early-early", "atom-atom", "se_atom-atom", "late-late", "se_late-late", "early-early minus late-late", "se_early-early minus late-late", "lo_early-early minus late-late", "hi_early-early minus late-late"]].round(3).to_string(index=False))
for k, v in pooled.items():
    print(k, {a: (round(b, 4) if isinstance(b, float) else b) for a, b in v.items() if a != "by_year"})
