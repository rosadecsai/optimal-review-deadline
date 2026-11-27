"""Figures 3 and 5 of the revised manuscript, redrawn from the archived outputs of the original pipeline
with the labels of the revised text (timing bins of Table 6; "final 48 h"; T* marked on V(T)).

usage: python code/revised_figures.py <repository root> <output folder>

Reads output/results.json. Nothing is re-estimated: Fig. 3 plots the archived cohort coefficients and Fig. 5
the delivery moments and values that code/revision_analyses.py reproduces exactly from the archived profile.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO, OUT = Path(sys.argv[1]), Path(sys.argv[2])
_argv = sys.argv; sys.argv = [_argv[0], str(REPO), str(OUT)]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import revision_analyses as RA            # definitions only (its analyses run under __main__)
sys.argv = _argv
RES = RA.RES
plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9, "legend.fontsize": 8, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5})
BIN_LABEL = {"mid": "2–4 days before", "last48h": "24–48 h before", "last24h": "final 24 h", "late": "late"}

# ---------------------------------------------------------------- Fig. 3: timing coefficients by cohort
by = pd.DataFrame(RES["E4"]["by_year"])
fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.4))
for ax, (out, ttl) in zip(axes, [("log_words", "(a) Length (log words)"), ("confidence", "(b) Confidence (1–5)"), ("dev_panel", "(c) Distance from co-reviewers' mean score")]):
    t = by[by["outcome"] == out]
    for (b, lab), m, dx in zip(BIN_LABEL.items(), ["o", "s", "^", "x"], [-0.18, -0.06, 0.06, 0.18]):
        tb = t[t["bin"] == b].dropna(subset=["coef"]).sort_values("year")
        ax.errorbar(tb["year"] + dx, tb["coef"], yerr=1.96 * tb["se"], marker=m, ms=3.5, lw=0.9, capsize=1.5, label=lab)
    ax.axhline(0, color="k", lw=0.7)
    ax.set_title(ttl); ax.set_xlabel("ICLR cohort"); ax.set_xticks(range(2018, 2027, 2))
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axes[0].set_ylabel("Coefficient relative to early reviews\n(four days or more before the deadline)")
axes[2].legend(frameon=False, loc="lower right")
fig.tight_layout(); fig.savefig(OUT / "fig3_quality_by_bin.pdf"); fig.savefig(OUT / "fig3_quality_by_bin.png", dpi=200); plt.close(fig)

# ---------------------------------------------------------------- Fig. 5: shares, delivery date and V(T)
sel = RA.TS <= 60; Ts = RA.TS[sel]
atom, late, _, _ = RA.PRE["power"]; atom, late = atom[sel], late[sel]
delay = RA.delay("power", RA.MU0)[0][sel]
val = pd.DataFrame(RES["E7"]["validation"])
fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.5))
ax = axes[0]
for y, lab, c, mk, col in ((atom, "final 48 h", "C0", "o", "obs_share_atom"), (late, "late", "C1", "s", "obs_share_late"), (1 - atom - late, "earlier", "C2", "^", "obs_share_pre")):
    ax.plot(Ts, 100 * y, color=c, label=lab)
    ax.plot(val["T"], 100 * val[col], mk, color=c, ms=3.5, mec="white", mew=0.4)
ax.set_xlabel("Review window T (days)"); ax.set_ylabel("Share of reviews (%)"); ax.set_title("(a) Predicted shares (markers: observed cohorts)")
ax.set_ylim(0, 100); ax.legend(frameon=False, ncol=3, loc="upper center")
ax = axes[1]
ax.plot(Ts, delay, color="C0", label="expected delivery")
ax.plot(Ts, Ts, ls="--", color="grey", lw=0.8, label="T")
ax.set_xlabel("Review window T (days)"); ax.set_ylabel("Days after assignment"); ax.set_title("(b) Expected delivery date"); ax.legend(frameon=False)
ax = axes[2]
aJ = (RA.ABAR0 * (1 - np.exp(-RA.LAM0 * RA.TS)))[None, :]; one = np.ones((1, len(RA.TS)))
for lab, acc, L, delta, eta in (("Conference, δ=0.005, η=0", one, 0.0, 0.005, 0.0), ("Conference, δ=0.005, η=1.3", one, 0.0, 0.005, 1.3),
                                ("Journal, L=0.3, δ=0.005, η=0", aJ, 0.3, 0.005, 0.0), ("Journal, L=0.3, δ=0.005, η=1.3", aJ, 0.3, 0.005, 1.3),
                                ("Journal, L=0.3, δ=0.01, η=0", aJ, 0.3, 0.01, 0.0)):
    V = RA.value(delta, L, RA.DA0, RA.DL0, eta, acc)[0][sel]; i = int(np.argmax(V))
    ln, = ax.plot(Ts, V, label=lab, lw=1.2)
    ax.plot([Ts[i]], [V[i]], "o", color=ln.get_color(), ms=4.5, mec="white", mew=0.6)
ax.set_xlabel("Review window T (days)"); ax.set_ylabel("Editor's expected value V(T)"); ax.set_title("(c) V(T) in five scenarios (dots: T*)")
ax.legend(frameon=False, fontsize=6.5, loc="lower center")
for a in axes:
    for s in ("top", "right"): a.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(OUT / "fig5_optimal_deadline.pdf"); fig.savefig(OUT / "fig5_optimal_deadline.png", dpi=200); plt.close(fig)
print("figures 3 and 5 written to", OUT)
