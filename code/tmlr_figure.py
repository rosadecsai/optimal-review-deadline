"""Figure: delivery relative to the deadline in the same-day panels of TMLR.
usage: python tmlr_figure.py <repository root> <tmlr output folder>"""
import csv, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO, OUT = sys.argv[1], sys.argv[2]
TM = json.load(open(os.path.join(OUT, "tmlr_results.json")))
RES = json.load(open(os.path.join(REPO, "output", "results.json")))
rows = list(csv.DictReader(open(os.path.join(OUT, "tmlr_same_day_panels_reviews.csv"))))
C14, C28, CREF = "#1f77b4", "#d95f02", "#666666"

plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "font.size": 9, "axes.labelsize": 9, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5})
fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.7))
ax = axes[0]
for g, col, lab, T in (("W14", C14, "14-day window", 14), ("W28", C28, "28-day window", 28)):
    t = np.sort([float(r["days_since_under_review"]) for r in rows if r["group"] == g])
    ax.step(t, np.arange(1, len(t) + 1) / len(t), where="post", color=col, lw=2)
    ax.axvline(T + 0.5, color=col, lw=0.9, ls=":")
    med = float(np.median(t))
    ax.plot([med], [0.5], "o", color=col, ms=6, mec="white", mew=1)
    if T == 14:
        ax.text(0.8, 0.56, f"{lab}\nmedian {med:.1f} d", fontsize=8, color="black", ha="left", va="bottom")
    else:
        ax.text(med + 1.6, 0.5, f"{lab}\nmedian {med:.1f} d", fontsize=8.5, color="black", ha="left", va="top")
ax.text(14.2, 1.02, "deadline", fontsize=8, color="black", ha="center"); ax.text(28.6, 1.02, "deadline", fontsize=8, color="black", ha="center")
ax.set_xlim(0, 60); ax.set_ylim(0, 1.0)
ax.set_xlabel("Days since the reviewers were assigned"); ax.set_ylabel("Share of reviews delivered")
ax.set_title("(a) Delivery by days since assignment", fontsize=10, pad=14)
ax.grid(axis="y", color="0.9", lw=0.6); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)

ax = axes[1]
for g, col, lab, mk in (("W14", C14, "TMLR, 14-day window", "o"), ("W28", C28, "TMLR, 28-day window", "s")):
    hz = {x["day"]: x for x in TM["hazard"][g]}
    xs, ys = [], []
    for a in (25, 22, 19, 16):                      # three-day bins beyond two weeks (few events per day)
        if all(d in hz for d in (a, a + 1, a + 2)):
            ev = sum(hz[d]["events"] for d in (a, a + 1, a + 2)); ar = sum(hz[d]["at_risk"] for d in (a, a + 1, a + 2))
            xs.append(a + 1); ys.append(ev / ar)
    for d in range(15 if g == "W28" else 13, -15, -1):
        if d in hz and hz[d]["hazard"] is not None and hz[d]["at_risk"] >= 30:
            xs.append(d); ys.append(hz[d]["hazard"])
    ax.plot(xs, ys, color=col, lw=1.8, marker=mk, ms=4, label=lab)
prof = RES["E7"]["profile"]; hi = {int(k): float(v) for k, v in prof["h"].items()}
ks = sorted(k for k in hi if k <= 14)
ax.plot(ks, [hi[k] for k in ks], color=CREF, lw=1.4, ls="--", label="ICLR 2024-2026 (conference)")
ax.plot([-1, -14], [float(prof["mu"])] * 2, color=CREF, lw=1.4, ls="--")
ax.axvline(-0.5, color="0.5", lw=0.8, ls=":"); ax.text(-0.9, 0.0043, "after the deadline", fontsize=8, ha="left", color="black")
ax.set_yscale("log"); ax.set_ylim(0.0035, 0.8); ax.set_xlim(27.5, -14.5)
ax.set_yticks([0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]); ax.set_yticklabels(["0.5%", "1%", "2%", "5%", "10%", "20%", "50%"])
ax.set_xlabel("Days before the deadline (negative: days late)"); ax.set_ylabel("Daily hazard of delivery (log scale)")
ax.set_title("(b) Hazard aligned on the deadline", fontsize=10, pad=14)
ax.legend(fontsize=8, frameon=False, loc="upper left"); ax.grid(axis="y", color="0.9", lw=0.6); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig4_tmlr_deadline.pdf")); fig.savefig(os.path.join(OUT, "fig4_tmlr_deadline.png"), dpi=170); plt.close(fig)
print("ok")
