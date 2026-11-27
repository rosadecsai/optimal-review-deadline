"""
Within-journal delivery relative to the review deadline: TMLR on OpenReview.

usage: python tmlr_analysis.py <tfd_tmlr_all.jsonl.gz> <tmlr_review_duedates.csv> <output folder>

Inputs
  tfd_tmlr_all.jsonl.gz      snapshot of the public TMLR submissions and their replies (13-14 Oct 2026),
                             from github.com/rosadecsai/time-to-first-decision (creation time of each review)
  tmlr_review_duedates.csv   due date of the per-paper review invitation TMLR/Paper<n>/-/Review
                             (read from the OpenReview API on 27 Nov 2026; papers put under review up to Dec 2025)

Design
  TMLR sets the review due date to the end (UTC) of the day 14 days after the reviewer assignment for regular
  submissions and 28 days after it for long submissions (more than 12 pages). The invitation keeps the date of the
  last assignment. A paper whose due date falls on the day 14 (28) days after the day it was put under review had
  all its reviewers assigned on that same day, so the whole panel shares one known deadline and one known window.
  These "same-day panels" are the analysis sample.
"""
import gzip, json, csv, sys, os, re, datetime
import numpy as np

DAY = 86400000.0
SEED = 20261127
NBOOT = 4000
BINS = ["early", "mid", "last48h", "last24h", "late"]

def tbin(s):
    return "late" if s < 0 else "last24h" if s < 1 else "last48h" if s < 2 else "mid" if s < 4 else "early"

def year(ms):
    return datetime.datetime.fromtimestamp(ms / 1000, datetime.UTC).year

def load(snap, duef):
    due = {int(r["number"]): (int(r["due_ms"]), int(r["end_of_day"])) for r in csv.DictReader(open(duef))}
    papers = []
    with gzip.open(snap, "rt", encoding="utf-8") as f:
        for line in f:
            d = json.loads(line)
            L = d.get("len") or ""
            cls = "R" if L.startswith("Regular") else "L" if L.startswith("Long") else "?"
            revs = sorted([x for x in d["rep"] if x["k"] == "R"], key=lambda x: x["c"])
            dec = [x for x in d["rep"] if x["k"] == "D"]
            papers.append(dict(num=d["num"], o=d["o"], cls=cls, v=d["v"], ae=d.get("ae"), revs=revs, dec=dec, due=due.get(d["num"])))
    return papers, due

def group_of(p, eps=1.0):
    """Window group of a paper, or None. eps = tolerated days between 'under review' and the last assignment.
    W14: due date 14 days after the day the paper was put under review (two-week window; regular at review time).
    W28: long submission whose due date falls 28 days after that day (four-week window).
    The length label in the snapshot is the current one: a paper reviewed as regular and relabelled long in a later
    revision has a 14-day window and is in W14 (label 'L').
    Papers with a review dated before the recorded under-review date are left out (the date was reset)."""
    if p["due"] is None or not p["o"] or p["cls"] not in "RL" or not p["revs"] or p["due"][1] != 1:
        return None
    if min(r["c"] for r in p["revs"]) < p["o"]:
        return None            # reviews dated before the recorded under-review date: the paper was put under review again
    w = (p["due"][0] - p["o"]) / DAY
    if 14 <= w < 14 + eps:
        return "W14"
    if p["cls"] == "L" and 28 <= w < 28 + eps:
        return "W28"
    return None

def reviews_of(papers, eps=1.0, years=None):
    rows = []
    for p in papers:
        g = group_of(p, eps)
        if g is None:
            continue
        if years and year(p["o"]) not in years:
            continue
        D = p["due"][0]
        for r in p["revs"]:
            rows.append(dict(num=p["num"], g=g, T=14 if g == "W14" else 28, lab=p["cls"], y=year(p["o"]), ae=p["ae"],
                             s=(D - r["c"]) / DAY, tau=(r["c"] - p["o"]) / DAY,
                             chars=r.get("len") or 0, ce=r.get("ce"), au=r.get("au"), nrev=len(p["revs"])))
    return rows

def q(a, p):
    return float(np.quantile(np.asarray(a, float), p))

def describe(rows):
    s = np.array([r["s"] for r in rows]); tau = np.array([r["tau"] for r in rows])
    n = len(rows); late = s[s < 0]
    out = dict(n_reviews=n, n_papers=len({r["num"] for r in rows}),
               tau_q25=q(tau, .25), tau_q50=q(tau, .5), tau_q75=q(tau, .75), tau_q90=q(tau, .9), tau_mean=float(tau.mean()),
               s_q10=q(s, .10), s_q25=q(s, .25), s_q50=q(s, .5), s_q75=q(s, .75), s_q90=q(s, .9),
               share_early=float((s >= 4).mean()), share_mid=float(((s >= 2) & (s < 4)).mean()),
               share_24_48h=float(((s >= 1) & (s < 2)).mean()), share_last24h=float(((s >= 0) & (s < 1)).mean()),
               share_final48h=float(((s >= 0) & (s < 2)).mean()), share_late=float((s < 0).mean()),
               share_first_half=float((s >= 14).mean()) if rows[0]["T"] == 28 else None,
               share_late_gt7=float((s < -7).mean()), share_late_gt14=float((s < -14).mean()), share_late_gt28=float((s < -28).mean()),
               late_median_days=float(np.median(-late)) if len(late) else None,
               ratio_24h_before_after=float(((s >= 0) & (s < 1)).sum() / max(1, ((s >= -1) & (s < 0)).sum())),
               ratio_48h_before_after=float(((s >= 0) & (s < 2)).sum() / max(1, ((s >= -2) & (s < 0)).sum())))
    return out

def hazard(rows, kmax_before, kmax_after=28):
    """daily hazard by day relative to the deadline: k = floor(s) (k>=0 before, k<0 after)"""
    s = np.array([r["s"] for r in rows]); k = np.floor(s).astype(int)
    out = []
    for d in list(range(kmax_before, -1, -1)) + list(range(-1, -kmax_after - 1, -1)):
        at_risk = int((k <= d).sum()); ev = int((k == d).sum())
        out.append(dict(day=d, at_risk=at_risk, events=ev, hazard=(ev / at_risk if at_risk else None)))
    return out

def cluster_index(rows):
    idx = {}
    for i, r in enumerate(rows):
        idx.setdefault(r["num"], []).append(i)
    return [np.array(v) for v in idx.values()]

def boot(rowsA, rowsB, stat, rng, nboot=NBOOT):
    """paper-clustered bootstrap of stat(valuesA, valuesB)"""
    cA, cB = cluster_index(rowsA), cluster_index(rowsB)
    out = np.empty(nboot)
    for b in range(nboot):
        ia = np.concatenate([cA[i] for i in rng.integers(0, len(cA), len(cA))])
        ib = np.concatenate([cB[i] for i in rng.integers(0, len(cB), len(cB))])
        out[b] = stat(ia, ib)
    return out

def passthrough(rows14, rows28, rng):
    t14 = np.array([r["tau"] for r in rows14]); t28 = np.array([r["tau"] for r in rows28])
    s14 = np.array([r["s"] for r in rows14]); s28 = np.array([r["s"] for r in rows28])
    res = {}
    for p in (.25, .5, .75, .9):
        f = lambda ia, ib, p=p: (np.quantile(t28[ib], p) - np.quantile(t14[ia], p)) / 14.0
        est = (np.quantile(t28, p) - np.quantile(t14, p)) / 14.0
        bs = boot(rows14, rows28, f, rng)
        res[f"q{int(p*100)}"] = dict(passthrough=float(est), lo=float(np.quantile(bs, .025)), hi=float(np.quantile(bs, .975)),
                                     slope_days_before=float(1 - est), tau14=float(np.quantile(t14, p)), tau28=float(np.quantile(t28, p)),
                                     s14=float(np.quantile(s14, 1 - p)), s28=float(np.quantile(s28, 1 - p)))
    f = lambda ia, ib: (t28[ib].mean() - t14[ia].mean()) / 14.0
    bs = boot(rows14, rows28, f, rng)
    res["mean"] = dict(passthrough=float((t28.mean() - t14.mean()) / 14.0), lo=float(np.quantile(bs, .025)), hi=float(np.quantile(bs, .975)),
                       tau14=float(t14.mean()), tau28=float(t28.mean()))
    # winsorised mean (delivery capped at 60 days after the deadline): robust to the few very late reports
    w14 = np.minimum(t14, 14 + 60); w28 = np.minimum(t28, 28 + 60)
    f = lambda ia, ib: (w28[ib].mean() - w14[ia].mean()) / 14.0
    bs = boot(rows14, rows28, f, rng)
    res["mean_capped60"] = dict(passthrough=float((w28.mean() - w14.mean()) / 14.0), lo=float(np.quantile(bs, .025)), hi=float(np.quantile(bs, .975)))
    return res

def share_diffs(rows14, rows28, rng):
    s14 = np.array([r["s"] for r in rows14]); s28 = np.array([r["s"] for r in rows28])
    defs = dict(final48h=lambda s: (s >= 0) & (s < 2), last24h=lambda s: (s >= 0) & (s < 1), late=lambda s: s < 0,
                early_ge4=lambda s: s >= 4, late_gt7=lambda s: s < -7)
    out = {}
    for k, fn in defs.items():
        a, b = fn(s14), fn(s28)
        bs = boot(rows14, rows28, lambda ia, ib: b[ib].mean() - a[ia].mean(), rng)
        out[k] = dict(share14=float(a.mean()), share28=float(b.mean()), diff=float(b.mean() - a.mean()),
                      lo=float(np.quantile(bs, .025)), hi=float(np.quantile(bs, .975)))
    return out

def common_support(rows14, rows28, rng, nperm=4000):
    """Proposition 1: conditional on not having delivered 14 days before the deadline, the 28-day reviewers should
    deliver like the 14-day reviewers. Kolmogorov-Smirnov distance on s (s<=14), paper-level permutation test."""
    a = [r for r in rows14 if r["s"] < 14]; b = [r for r in rows28 if r["s"] < 14]
    def ks(x, y):
        x = np.sort(x); y = np.sort(y); allv = np.concatenate([x, y])
        return float(np.max(np.abs(np.searchsorted(x, allv, side="right") / len(x) - np.searchsorted(y, allv, side="right") / len(y))))
    sa = np.array([r["s"] for r in a]); sb = np.array([r["s"] for r in b])
    obs = ks(sa, sb)
    # permutation of the group label at paper level
    allr = a + b; cl = cluster_index(allr); s_all = np.array([r["s"] for r in allr])
    nA = len({r["num"] for r in a}); cnt = 0
    for _ in range(nperm):
        perm = rng.permutation(len(cl)); ia = np.concatenate([cl[i] for i in perm[:nA]]); ib = np.concatenate([cl[i] for i in perm[nA:]])
        cnt += ks(s_all[ia], s_all[ib]) >= obs
    d = lambda s: dict(n=int(len(s)), final48h=float(((s >= 0) & (s < 2)).mean()), late=float((s < 0).mean()), median=float(np.median(s)),
                       q25=float(np.quantile(s, .25)), q75=float(np.quantile(s, .75)))
    return dict(ks=obs, p_perm=float((cnt + 1) / (nperm + 1)), cond14=d(sa), cond28=d(sb),
                share28_delivered_before_last14=float(np.mean([r["s"] >= 14 for r in rows28])))

def fe_reg(rows, yname, transform=None):
    """paper fixed effects, timing bins (reference: early), paper-clustered SE"""
    rows = [r for r in rows if r.get(yname) is not None]
    if transform:
        rows = [dict(r, **{yname: transform(r[yname])}) for r in rows if transform(r[yname]) is not None]
    cl = [c for c in cluster_index(rows) if len(c) >= 2]
    keep = np.concatenate(cl); rows = [rows[i] for i in keep]; cl = cluster_index(rows)
    y = np.array([r[yname] for r in rows], float)
    X = np.array([[tbin(r["s"]) == b for b in BINS[1:]] for r in rows], float)
    for c in cl:
        y[c] -= y[c].mean(); X[c] -= X[c].mean(axis=0)
    XtXi = np.linalg.pinv(X.T @ X); beta = XtXi @ X.T @ y; e = y - X @ beta
    meat = np.zeros((X.shape[1], X.shape[1]))
    for c in cl:
        g = X[c].T @ e[c]; meat += np.outer(g, g)
    G = len(cl); n, k = X.shape
    V = XtXi @ meat @ XtXi * (G / (G - 1)) * ((n - 1) / (n - k - G))
    se = np.sqrt(np.diag(V))
    return dict(n=int(n), papers=int(G), coef={b: dict(b=float(beta[i]), se=float(se[i])) for i, b in enumerate(BINS[1:])},
                bin_n={b: int(sum(tbin(r["s"]) == b for r in rows)) for b in BINS})

def main(snap, duef, outdir):
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(SEED)
    papers, due = load(snap, duef)
    R = {}
    # ---------- sample accounting
    with_due = [p for p in papers if p["due"] is not None and p["o"] and p["cls"] in "RL"]
    acc = dict(snapshot_submissions=len(papers), due_dates_read=len(due), with_due_and_class=len(with_due),
               with_due_and_reviews=sum(1 for p in with_due if p["revs"]),
               regular=sum(p["cls"] == "R" for p in with_due), long=sum(p["cls"] == "L" for p in with_due),
               reviews_all=sum(len(p["revs"]) for p in with_due))
    w = {c: np.array([(p["due"][0] - p["o"]) / DAY for p in with_due if p["cls"] == c]) for c in "RL"}
    acc["regular_window_min"] = float(w["R"][w["R"] > 0].min()); acc["regular_share_lt14"] = float((w["R"] < 14).mean())
    acc["regular_share_same_day"] = float(((w["R"] >= 14) & (w["R"] < 15)).mean())
    acc["long_share_same_day_28"] = float(((w["L"] >= 28) & (w["L"] < 29)).mean())
    acc["long_share_same_day_14"] = float(((w["L"] >= 14) & (w["L"] < 15)).mean())
    acc["long_share_lt14"] = float((w["L"] < 14).mean())
    acc["median_days_under_review_to_due_regular"] = float(np.median(w["R"])); acc["median_days_under_review_to_due_long"] = float(np.median(w["L"]))
    # reviews that precede the last assignment (due - window - 1 day), by distance between 'under review' and the due date
    pre = {}
    for lo, hi in ((-1e9, 0), (0, 1), (1, 7), (7, 21), (21, 60), (60, 1e9)):
        n = k = npap = 0
        for p in with_due:
            T = 14 if p["cls"] == "R" else 28; sl = (p["due"][0] - p["o"]) / DAY - T
            if lo <= sl < hi and p["revs"]:
                npap += 1
                for r in p["revs"]:
                    n += 1; k += (p["due"][0] - r["c"]) / DAY > T + 1
        key = "neg" if lo < 0 else f"{lo}-{hi if hi < 1e8 else 'inf'}"
        pre[key] = dict(papers=npap, reviews=n, before_last_assignment=k, share=(k / n if n else None))
    acc["reviews_before_last_assignment_by_slack"] = pre
    acc["with_reviews_due_date_not_end_of_day"] = sum(1 for p in with_due if p["revs"] and p["due"][1] != 1)
    def _raw_group(p):
        if p["due"] is None or not p["o"] or p["cls"] not in "RL" or not p["revs"] or p["due"][1] != 1:
            return None
        w_ = (p["due"][0] - p["o"]) / DAY
        return "W14" if 14 <= w_ < 15 else "W28" if (p["cls"] == "L" and 28 <= w_ < 29) else None
    dropped = [p for p in with_due if _raw_group(p) and min(r["c"] for r in p["revs"]) < p["o"]]
    acc["dropped_reviews_before_under_review"] = dict(papers=len(dropped), reviews=sum(len(p["revs"]) for p in dropped),
                                                     days_before=[round(float((p["o"] - min(r["c"] for r in p["revs"])) / DAY), 1) for p in dropped])
    R["accounting"] = acc
    # ---------- main sample: same-day panels, papers put under review 2022-2025
    rows = reviews_of(papers, 1.0)
    G = {g: [r for r in rows if r["g"] == g] for g in ("W14", "W28")}
    G["W14_regular_label"] = [r for r in G["W14"] if r["lab"] == "R"]
    G["W14_long_label"] = [r for r in G["W14"] if r["lab"] == "L"]
    R["groups"] = {g: describe(v) for g, v in G.items()}
    R["n_reviews_per_paper"] = {}
    for g in ("W14", "W28"):
        cnt = {}
        for c in cluster_index(G[g]):
            cnt[len(c)] = cnt.get(len(c), 0) + 1
        R["n_reviews_per_paper"][g] = {str(k): v for k, v in sorted(cnt.items())}
    R["groups_by_year"] = {g: {str(y): dict(n_papers=len({r["num"] for r in G[g] if r["y"] == y}), n_reviews=sum(r["y"] == y for r in G[g])) for y in sorted({r["y"] for r in G[g]})} for g in ("W14", "W28")}
    R["hazard"] = {"W14": hazard(G["W14"], 13), "W28": hazard(G["W28"], 27)}
    # pooled pre-deadline profile on the common support (used by the calibration)
    pooled = G["W14"] + [r for r in G["W28"] if r["s"] < 14]
    R["hazard"]["pooled_last14"] = hazard(pooled, 13)
    R["passthrough"] = passthrough(G["W14"], G["W28"], rng)
    R["share_diffs"] = share_diffs(G["W14"], G["W28"], rng)
    R["common_support"] = common_support(G["W14"], G["W28"], rng)
    R["relabelled_vs_regular_within_W14"] = dict(common_support(G["W14_regular_label"], G["W14_long_label"], rng),
                                                 passthrough_like_q50=float((np.median([r["tau"] for r in G["W14_long_label"]]) - np.median([r["tau"] for r in G["W14_regular_label"]]))))
    def tail(rows_, lo, hi):
        s = np.array([r["s"] for r in rows_]); k = np.floor(s).astype(int); ev = ar = 0
        for d in range(-lo, -hi - 1, -1):
            ar += int((k <= d).sum()); ev += int((k == d).sum())
        return dict(events=ev, at_risk_days=ar, hazard=(ev / ar if ar else None))
    R["late_hazard"] = {g: {f"{lo}-{hi}": tail(G[g], lo, hi) for lo, hi in ((1, 7), (8, 14), (15, 28), (1, 28), (1, 60))} for g in ("W14", "W28")}
    R["late_hazard"]["pooled"] = {f"{lo}-{hi}": tail(G["W14"] + G["W28"], lo, hi) for lo, hi in ((1, 7), (8, 14), (15, 28), (1, 28), (1, 60))}
    # pooled daily hazard before the deadline (days 13..0) and in the first week after it: difference between windows
    def pooled_h(sv, lo, hi):
        k = np.floor(sv).astype(int); ev = ar = 0
        for d in range(hi, lo - 1, -1):
            ar += int((k <= d).sum()); ev += int((k == d).sum())
        return ev / ar if ar else np.nan
    s14 = np.array([r["s"] for r in G["W14"]]); s28 = np.array([r["s"] for r in G["W28"]])
    hc = {}
    for name, lo, hi in (("before_13_0", 0, 13), ("after_1_7", -7, -1), ("after_1_14", -14, -1)):
        bs = boot(G["W14"], G["W28"], lambda ia, ib, lo=lo, hi=hi: pooled_h(s28[ib], lo, hi) - pooled_h(s14[ia], lo, hi), rng, 2000)
        hc[name] = dict(h14=float(pooled_h(s14, lo, hi)), h28=float(pooled_h(s28, lo, hi)), diff=float(pooled_h(s28, lo, hi) - pooled_h(s14, lo, hi)),
                        lo=float(np.quantile(bs, .025)), hi=float(np.quantile(bs, .975)))
    R["hazard_contrast"] = hc
    # action editors: bootstrap that resamples action editors (with all their panels), and the comparison within editors
    allr = G["W14"] + G["W28"]; tau_all = np.array([r["tau"] for r in allr]); is28 = np.array([r["g"] == "W28" for r in allr])
    by_ae = {}
    for i, r in enumerate(allr):
        by_ae.setdefault(r["ae"], []).append(i)
    cl = [np.array(v) for v in by_ae.values()]; bs = []
    for _ in range(NBOOT):
        ii = np.concatenate([cl[j] for j in rng.integers(0, len(cl), len(cl))])
        a = tau_all[ii][~is28[ii]]; b = tau_all[ii][is28[ii]]
        if len(a) and len(b): bs.append((np.median(b) - np.median(a)) / 14.0)
    ae14 = {r["ae"] for r in G["W14"]}; ae28 = {r["ae"] for r in G["W28"]}; both = ae14 & ae28
    dm = []; wt = []
    for a_ in both:
        ta = [r["tau"] for r in G["W14"] if r["ae"] == a_]; tb = [r["tau"] for r in G["W28"] if r["ae"] == a_]
        dm.append(np.median(tb) - np.median(ta)); wt.append(len(ta) * len(tb) / (len(ta) + len(tb)))
    R["action_editors"] = dict(n_editors=len(by_ae), n_editors_both_windows=len(both),
                               passthrough_q50_editor_bootstrap=dict(lo=float(np.quantile(bs, .025)), hi=float(np.quantile(bs, .975))),
                               within_editor_weighted_median_diff_over_14=float(np.average(dm, weights=wt) / 14.0),
                               within_editor_median_of_diffs_over_14=float(np.median(dm) / 14.0),
                               papers_in_both=dict(W14=len({r["num"] for r in G["W14"] if r["ae"] in both}), W28=len({r["num"] for r in G["W28"] if r["ae"] in both})))
    # elapsed-time hazards on the day the short window ends: 14-day reviewers vs 28-day reviewers at the same elapsed time
    def elapsed_hazard(rows_, d0, d1):
        t = np.array([r["tau"] for r in rows_]); ev = int(((t >= d0) & (t < d1)).sum()); ar = int((t >= d0).sum())
        return dict(events=ev, at_risk=ar, share=(ev / ar if ar else None))
    R["elapsed_days_12_15"] = {g: elapsed_hazard(G[g], 12.0, 15.0) for g in ("W14", "W28")}
    R["elapsed_days_26_29"] = {g: elapsed_hazard(G[g], 26.0, 29.0) for g in ("W14", "W28")}
    # ---------- report length and the claims-and-evidence answer by timing (paper fixed effects)
    logc = lambda c: float(np.log(c)) if c and c > 0 else None
    allrows = G["W14"] + G["W28"]
    R["fe_log_chars"] = {"all": fe_reg(allrows, "chars", logc), "W14": fe_reg(G["W14"], "chars", logc), "W28": fe_reg(G["W28"], "chars", logc)}
    def fe_single(rows_, yname, transform, flag):
        """paper fixed effects with one indicator; paper-clustered SE"""
        rows_ = [dict(r, y=transform(r[yname])) for r in rows_ if r.get(yname) is not None and transform(r[yname]) is not None]
        cl = [c for c in cluster_index(rows_) if len(c) >= 2]; rows_ = [rows_[i] for i in np.concatenate(cl)]; cl = cluster_index(rows_)
        y = np.array([r["y"] for r in rows_], float); x = np.array([flag(r) for r in rows_], float)
        for c in cl:
            y[c] -= y[c].mean(); x[c] -= x[c].mean()
        b = float((x @ y) / (x @ x)); e = y - b * x
        v = sum(float(x[c] @ e[c]) ** 2 for c in cl) / (x @ x) ** 2 * (len(cl) / (len(cl) - 1))
        return dict(b=b, se=float(np.sqrt(v)), n=len(rows_), papers=len(cl), share_flag=float(np.mean([flag(r) for r in rows_])))
    R["fe_log_chars_late_vs_on_time"] = fe_single(allrows, "chars", logc, lambda r: r["s"] < 0)
    R["fe_log_chars_final48h_vs_earlier_on_time"] = fe_single([r for r in allrows if r["s"] >= 0], "chars", logc, lambda r: r["s"] < 2)
    R["fe_log_chars_per_week_late"] = fe_single(allrows, "chars", logc, lambda r: max(0.0, -r["s"]) / 7.0)
    yes = lambda v: None if v is None else 1.0 if str(v).startswith("Yes") else 0.0
    R["fe_claims_supported"] = {"all": fe_reg([r for r in allrows if r["ce"] is not None], "ce", yes)}
    R["fe_claims_supported_late_vs_on_time"] = fe_single([r for r in allrows if r["ce"] is not None], "ce", yes, lambda r: r["s"] < 0)
    R["claims_supported_share"] = float(np.mean([str(r["ce"]).startswith("Yes") for r in allrows if r["ce"] is not None]))
    R["chars_median"] = {g: float(np.median([r["chars"] for r in G[g] if r["chars"]])) for g in ("W14", "W28")}
    # between-window difference in report length (papers differ in length: descriptive)
    l14 = np.array([np.log(r["chars"]) if r["chars"] else np.nan for r in G["W14"]]); l28 = np.array([np.log(r["chars"]) if r["chars"] else np.nan for r in G["W28"]])
    bs = boot(G["W14"], G["W28"], lambda ia, ib: np.nanmean(l28[ib]) - np.nanmean(l14[ia]), rng)
    R["log_chars_W28_minus_W14"] = dict(diff=float(np.nanmean(l28) - np.nanmean(l14)), lo=float(np.quantile(bs, .025)), hi=float(np.quantile(bs, .975)))
    # ---------- lower bound of lateness in the full sample: reviews submitted after the final due date of the paper
    fb = {}
    for c in ("R", "L", "RL"):
        n = k = k7 = 0
        for p in with_due:
            if p["cls"] in c:
                for r in p["revs"]:
                    n += 1; k += r["c"] > p["due"][0]; k7 += r["c"] > p["due"][0] + 7 * DAY
        fb[c] = dict(reviews=n, after_final_due=k, share=k / n, share_gt7=k7 / n)
    R["full_sample_after_final_due"] = fb
    # ---------- robustness
    rob = {}
    for name, eps, yrs in (("eps2", 2.0, None), ("eps3", 3.0, None), ("y2023_2025", 1.0, {2023, 2024, 2025}), ("y2024_2025", 1.0, {2024, 2025}), ("y2022_2023", 1.0, {2022, 2023})):
        rr = reviews_of(papers, eps, yrs); a = [r for r in rr if r["g"] == "W14"]; b = [r for r in rr if r["g"] == "W28"]
        pt = passthrough(a, b, rng)
        rob[name] = dict(W14=describe(a), W28=describe(b), passthrough_q50=pt["q50"], passthrough_q25=pt["q25"], passthrough_q75=pt["q75"], passthrough_mean=pt["mean"])
    three = [r for r in rows if r["nrev"] == 3]
    a = [r for r in three if r["g"] == "W14"]; b = [r for r in three if r["g"] == "W28"]
    pt = passthrough(a, b, rng)
    rob["exactly_three_reviews"] = dict(W14=describe(a), W28=describe(b), passthrough_q50=pt["q50"], passthrough_q25=pt["q25"], passthrough_q75=pt["q75"], passthrough_mean=pt["mean"])
    a = G["W14_regular_label"]; pt = passthrough(a, G["W28"], rng)
    rob["W14_regular_label_only"] = dict(W14=describe(a), W28=describe(G["W28"]), passthrough_q50=pt["q50"], passthrough_q25=pt["q25"], passthrough_q75=pt["q75"], passthrough_mean=pt["mean"])
    R["robustness"] = rob
    # ---------- how selected are the same-day panels?
    def paper_stats(ps):
        first = [(p["revs"][0]["c"] - p["o"]) / DAY for p in ps]; last = [(p["revs"][-1]["c"] - p["o"]) / DAY for p in ps]
        acc_ = [any(str(d.get("rec") or "").startswith("Accept") for d in p["dec"]) for p in ps if p["dec"]]
        return dict(n=len(ps), n_reviews=sum(len(p["revs"]) for p in ps), decided=len(acc_),
                    first_review_median=float(np.median(first)), last_review_median=float(np.median(last)),
                    mean_reviews=float(np.mean([len(p["revs"]) for p in ps])), accept_share_of_decided=float(np.mean(acc_)) if acc_ else None,
                    chars_median=float(np.median([r.get("len") or 0 for p in ps for r in p["revs"]])))
    sel = {}
    for c in "RL":
        base = [p for p in with_due if p["cls"] == c and p["revs"]]
        same = [p for p in base if group_of(p) == ("W14" if c == "R" else "W28")]
        other = [p for p in base if group_of(p) is None]
        sel[c] = dict(same_day=paper_stats(same), other=paper_stats(other))
    for c in "RL":
        base = [p for p in with_due if p["cls"] == c and p["revs"]]
        sel[c]["share_same_day"] = sel[c]["same_day"]["n"] / len(base); sel[c]["n_all_with_reviews"] = len(base)
        a_, b_ = sel[c]["same_day"], sel[c]["other"]; pa, pb = a_["accept_share_of_decided"], b_["accept_share_of_decided"]
        se_ = (pa * (1 - pa) / a_["decided"] + pb * (1 - pb) / b_["decided"]) ** 0.5
        sel[c]["accept_diff"] = dict(diff=pa - pb, se=se_, z=(pa - pb) / se_)
        yrs = sorted({year(p["o"]) for p in base})
        sel[c]["share_same_day_by_year"] = {str(y): (sum(1 for p in base if year(p["o"]) == y and group_of(p) == ("W14" if c == "R" else "W28")) / max(1, sum(1 for p in base if year(p["o"]) == y))) for y in yrs}
    R["selection"] = sel
    R["meta"] = dict(seed=SEED, nboot=NBOOT, snapshot="tfd_tmlr_all.jsonl.gz (13-14 Oct 2026)", due_dates="read 27 Nov 2026")
    json.dump(R, open(os.path.join(outdir, "tmlr_results.json"), "w"), indent=1)
    # per-review file of the analysis sample (no text, no identities)
    with open(os.path.join(outdir, "tmlr_same_day_panels_reviews.csv"), "w", newline="") as f:
        wcsv = csv.writer(f); wcsv.writerow(["paper", "group", "window_days", "year_under_review", "days_before_deadline", "days_since_under_review", "review_chars"])
        for r in rows:
            wcsv.writerow([r["num"], r["g"], r["T"], r["y"], f"{r['s']:.4f}", f"{r['tau']:.4f}", r["chars"]])
    with open(os.path.join(outdir, "tmlr_hazard_by_day.csv"), "w", newline="") as f:
        wcsv = csv.writer(f); wcsv.writerow(["group", "day_relative_to_deadline", "at_risk", "events", "hazard"])
        for g, h in R["hazard"].items():
            for x in h:
                wcsv.writerow([g, x["day"], x["at_risk"], x["events"], "" if x["hazard"] is None else f"{x['hazard']:.5f}"])
    return R

if __name__ == "__main__":
    R = main(sys.argv[1], sys.argv[2], sys.argv[3])
    rd = lambda o: json.dumps(o, default=float)
    for g, d in R["groups"].items():
        print(g, json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in d.items()}))
    for k in ("hazard_contrast", "action_editors", "n_reviews_per_paper", "groups_by_year", "passthrough", "share_diffs", "common_support", "relabelled_vs_regular_within_W14", "late_hazard", "elapsed_days_12_15", "elapsed_days_26_29", "fe_log_chars", "fe_log_chars_late_vs_on_time", "fe_log_chars_final48h_vs_earlier_on_time", "fe_log_chars_per_week_late", "fe_claims_supported_late_vs_on_time", "claims_supported_share", "chars_median", "log_chars_W28_minus_W14", "full_sample_after_final_due", "selection"):
        print(k, re.sub(r"(\d+\.\d{4})\d+", r"\1", rd(R[k]))[:2600])
    for k, v in R["robustness"].items():
        print(k, v["W14"]["n_papers"], v["W28"]["n_papers"], "q50", {x: round(v["passthrough_q50"][x], 3) for x in ("passthrough", "lo", "hi")}, "mean", {x: round(v["passthrough_mean"][x], 3) for x in ("passthrough", "lo", "hi")}, "48h", round(v["W14"]["share_final48h"], 3), round(v["W28"]["share_final48h"], 3), "late", round(v["W14"]["share_late"], 3), round(v["W28"]["share_late"], 3))
