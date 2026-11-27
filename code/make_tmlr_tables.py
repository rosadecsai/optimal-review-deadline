"""LaTeX tables of the TMLR analysis.
usage: python make_tmlr_tables.py <tmlr output folder> <manuscript tables folder> <first supplementary table number>
Writes tables/tab_tmlr.tex (main text) and, next to the tables folder, OnlineResource3_tmlr.tex (section of Online Resource 3)."""
import json, sys, os, csv, math
OUT, TAB, N0 = sys.argv[1], sys.argv[2], int(sys.argv[3])
R = json.load(open(os.path.join(OUT, "tmlr_results.json"))); C = json.load(open(os.path.join(OUT, "tmlr_calibration.json")))
G = R["groups"]; PT = R["passthrough"]
pc = lambda x, d=0: f"{100*x:.{d}f}"
f2 = lambda x: f"{x:.2f}"; f1 = lambda x: f"{x:.1f}"
num = lambda n: f"{n:,}".replace(",", "{,}")
sgn = lambda x, d=2: (f"${x:.{d}f}$" if x >= 0 else f"$-{abs(x):.{d}f}$")
ci = lambda d, k=2: f"{d['passthrough']:.{k}f} ({d['lo']:.{k}f}--{d['hi']:.{k}f})"

# ------------------------------------------------------------------ main-text table
rows = []
for g, lab in (("W14", "14 days"), ("W28", "28 days")):
    d = G[g]
    rows.append(f"{lab} & {num(d['n_papers'])} & {num(d['n_reviews'])} & {f1(d['tau_q50'])} & {sgn(d['s_q50'],1)} & {pc(1-d['share_final48h']-d['share_late'])} & {pc(d['share_final48h'])} & {pc(d['share_last24h'])} & {pc(d['share_late'])} & {f1(d['late_median_days'])} & {pc(d['share_late_gt7']/d['share_late'])} \\\\")
main = r"""\begin{table}[t]
\caption{\rev{One journal, two windows: delivery relative to the deadline in the same-day panels of TMLR}\revonly{ (new table)}}\label{tab:tmlr}
\scriptsize
\setlength{\tabcolsep}{2pt}
\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lrrrrrrrrrr@{}}
\toprule
Window & Papers & Reviews & \multicolumn{2}{c}{Median days} & \multicolumn{4}{c}{Share of the reviews (\%)} & \multicolumn{2}{c}{Late reviews} \\
 & & & since & before & $>$48 h & final & (final & late & median & $>$7 d \\
 & & & assignment & deadline & before & 48 h & 24 h) & & days late & (\%) \\
\midrule
""" + "\n".join(rows) + r"""
\midrule
\multicolumn{11}{@{}p{\textwidth}@{}}{\emph{Pass-through of the window to the delivery date} (difference between the two windows in days since assignment, divided by 14; 95\% interval from """ + num(R["meta"]["nboot"]) + r""" bootstrap resamples of papers): lower quartile """ + ci(PT["q25"]) + r"""; median """ + ci(PT["q50"]) + r"""; upper quartile """ + ci(PT["q75"]) + r"""; ninth decile """ + ci(PT["q90"]) + r"""; mean """ + ci(PT["mean"]) + r""". The intervals measure sampling among papers and take the two groups of manuscripts as given.} \\
\bottomrule
\end{tabular*}
\footnotetext{Transactions on Machine Learning Research. The 14-day window applies to regular submissions, of at most twelve pages of main content, and the 28-day window to long ones. Papers put under review from 2022 to 2025 whose reviewers were all assigned on the day the paper was put under review, so that every reviewer of the panel has the same window and the same deadline (the end of the day 14 or 28 days later, UTC). The window is set by the length of the manuscript, not at random. The 14-day panels include """ + str(G['W14_long_label']['n_papers']) + r""" submissions that carry the label long today; in those whose history we read, the label was changed in a revision after the reviews (Online Resource~3). Columns as in Table~\ref{tab:delivery}: the shares more than 48 hours before, in the final 48 hours and late sum to 100; a negative median of days before the deadline is a median after it; lateness statistics over late reviews, and the last column is a share of late reviews. A pass-through of 1 means that delivery moves one day later for each day added to the window.}
\end{table}
"""
open(os.path.join(TAB, "tab_tmlr.tex"), "w").write(main)

# ------------------------------------------------------------------ Online Resource 3, TMLR section
n = N0; S = []
def table(caption, label, colspec, header, body, note, size="\\footnotesize", env="table", width=None):
    """article-class table of Online Resource 3 (note in a minipage below the table)"""
    global n
    mp = "0.95\\textheight" if env == "sidewaystable" else "0.95\\textwidth"
    s = (f"\\begin{{{env}}}[htbp]\n\\caption{{{caption}}}\\label{{{label}}}\n\\centering{size}\\setlength{{\\tabcolsep}}{{4pt}}\n\\begin{{tabular}}{{{colspec}}}\n\\toprule\n{header} \\\\\n\\midrule\n"
         + "\n".join(body) + f"\n\\bottomrule\n\\end{{tabular}}\n\n\\smallskip\n\\begin{{minipage}}{{{mp}}}\\footnotesize {note}\\end{{minipage}}\n\\end{{{env}}}\n")
    n += 1; return s

A = R["accounting"]; pre = A["reviews_before_last_assignment_by_slack"]
T_acc = N0
body = [f"Public submissions in the snapshot (13--14 October 2026) & {num(A['snapshot_submissions'])} \\\\",
        f"Review due dates read from OpenReview (papers put under review up to December 2025) & {num(A['due_dates_read'])} \\\\",
        f"\\quad with a length label and an under-review date in the snapshot & {num(A['with_due_and_class'])} \\\\",
        f"\\quad\\quad of which with public reviews & {num(A['with_due_and_reviews'])} ({num(A['reviews_all'])} reviews) \\\\",
        f"\\quad\\quad of which with a due date not at the end of a day (set by hand; left out) & {num(A['with_reviews_due_date_not_end_of_day'])} \\\\",
        f"Panels with all reviews dated before the recorded under-review date (date reset; left out) & {A['dropped_reviews_before_under_review']['papers']} papers, {A['dropped_reviews_before_under_review']['reviews']} reviews \\\\",
        f"Same-day panels, 14-day window & {num(G['W14']['n_papers'])} papers, {num(G['W14']['n_reviews'])} reviews \\\\",
        f"\\quad label ``regular'' / label ``long'' (relabelled in a later revision) & {num(G['W14_regular_label']['n_papers'])} / {num(G['W14_long_label']['n_papers'])} papers \\\\",
        f"Same-day panels, 28-day window & {num(G['W28']['n_papers'])} papers, {num(G['W28']['n_reviews'])} reviews \\\\",
        "\\midrule",
        "\\multicolumn{2}{@{}p{0.9\\textwidth}@{}}{\\emph{All " + num(A['reviews_all']) + " reviews: those submitted before the last assignment (more than window $+1$ days before the final due date), by days between the under-review date and the last assignment; window by the current length label}} \\\\"]
for k, lab in (("0-1", "same day"), ("1-7", "1--7 days"), ("7-21", "7--21 days"), ("21-60", "21--60 days"), ("60-inf", "more than 60 days"), ("neg", "due date earlier than the current label allows (relabelled after review)")):
    body.append(f"\\quad {lab} & {num(pre[k]['before_last_assignment'])} of {num(pre[k]['reviews'])} ({pc(pre[k]['share'],1)}\\%) \\\\")
fb = R["full_sample_after_final_due"]["RL"]
body += ["\\midrule", f"Reviews submitted after the final due date of their paper, all {num(fb['reviews'])} reviews & {pc(fb['share'],1)}\\% (more than 7 days after: {pc(fb['share_gt7'],1)}\\%) \\\\"]
S.append(table("TMLR: construction of the same-day panels", "tab:tmlr-acc", "p{0.62\\textwidth}r", "Step & Count", body,
               "The workflow of the journal sets the due date of the review invitation of a paper to the end of the day (UTC) two weeks (regular submissions) or four weeks (long submissions) after the paper is put under review, and resets it in the same way each time a reviewer is assigned; the invitation therefore keeps the date of the last assignment. A paper whose due date falls on the day 14 (28) days after the day it was put under review had no reviewer assigned after that day. The increasing share of reviews that precede the last assignment (middle block) confirms this reading: such reviews can only come from reviewers assigned earlier. That block classifies every paper by its current length label, so its first row is not the analysis sample: it contains the three panels left out above and a few papers with a due date set by hand, and the panels relabelled long after the review are in its last row. The last row of the table is a lower bound of lateness over all reviews, because the final due date is the latest deadline any reviewer of the paper had."))

T_del = n
hdr = "Group & Papers & Reviews & \\multicolumn{3}{c}{Days since assignment} & \\multicolumn{3}{c}{Days before the deadline} & $\\ge$4 d & 2--4 d & 24--48 h & Final 24 h & Late & Late $>$7 d & $>$14 d & $>$28 d \\\\\n & & & Q25 & Q50 & Q75 & Q25 & Q50 & Q75 & \\multicolumn{8}{c}{(\\% of the reviews of the group)}"
body = []
for g, lab in (("W14", "14-day window"), ("W14_regular_label", "\\quad label regular"), ("W14_long_label", "\\quad label long"), ("W28", "28-day window")):
    d = G[g]
    body.append(f"{lab} & {num(d['n_papers'])} & {num(d['n_reviews'])} & {f1(d['tau_q25'])} & {f1(d['tau_q50'])} & {f1(d['tau_q75'])} & {sgn(d['s_q25'],1)} & {sgn(d['s_q50'],1)} & {sgn(d['s_q75'],1)} & {pc(d['share_early'],1)} & {pc(d['share_mid'],1)} & {pc(d['share_24_48h'],1)} & {pc(d['share_last24h'],1)} & {pc(d['share_late'],1)} & {pc(d['share_late_gt7'],1)} & {pc(d['share_late_gt14'],1)} & {pc(d['share_late_gt28'],1)} \\\\")
cs = R["common_support"]; rl = R["relabelled_vs_regular_within_W14"]
S.append(table("TMLR: delivery by window in the same-day panels", "tab:tmlr-del", "lrrrrrrrrrrrrrrrr", hdr, body,
               f"Quantiles of days before the deadline are those of the distribution of $s$ (negative: after the deadline), so the lower quartile is the late side. All shares over the reviews of the group. Of the reviews with a 28-day window, {pc(cs['share28_delivered_before_last14'],1)}\\% were delivered in the first 14 days. Restricted to the last 14 days before the deadline and afterwards, the two windows differ (largest distance between the two cumulative distributions, the Kolmogorov--Smirnov statistic, {cs['ks']:.3f}; permutation of the window label over papers $p={cs['p_perm']:.4f}$; {num(cs['cond14']['n'])} and {num(cs['cond28']['n'])} reviews): the 28-day reviewers who reach the last two weeks are late somewhat more often ({pc(cs['cond28']['late'],1)}\\% against {pc(cs['cond14']['late'],1)}\\%). Within the 14-day window, papers labelled regular and papers relabelled long do not differ (distance {rl['ks']:.3f}, $p={rl['p_perm']:.2f}$).",
               size="\\scriptsize", env="sidewaystable"))

T_haz = n
hz = {g: {x["day"]: x for x in R["hazard"][g]} for g in ("W14", "W28")}
body = []
for d in list(range(27, -1, -1)) + list(range(-1, -15, -1)):
    a = hz["W14"].get(d); b = hz["W28"].get(d)
    fa = f"{num(a['at_risk'])} & {a['events']} & {pc(a['hazard'],1)}" if a else " & & "
    fb_ = f"{num(b['at_risk'])} & {b['events']} & {pc(b['hazard'],1)}" if b else " & & "
    body.append(f"{d if d>=0 else '$-'+str(-d)+'$'} & {fa} & {fb_} \\\\")
lh = R["late_hazard"]["pooled"]; l14 = R["late_hazard"]["W14"]; l28 = R["late_hazard"]["W28"]; hc = R["hazard_contrast"]
S.append(table("TMLR: daily hazard of delivery by day relative to the deadline", "tab:tmlr-haz", "rrrrrrr",
               "Days before the deadline & \\multicolumn{3}{c}{14-day window} & \\multicolumn{3}{c}{28-day window} \\\\\n & At risk & Delivered & Hazard (\\%) & At risk & Delivered & Hazard (\\%)", body,
               f"Day 0 is the final 24 hours before the deadline; negative days are days after it. At risk: reviews of the group not yet delivered at the start of the day. Pooled over the two windows, the daily hazard after the deadline is {lh['1-7']['hazard']:.3f} in days 1--7, {lh['8-14']['hazard']:.3f} in days 8--14, {lh['15-28']['hazard']:.3f} in days 15--28 and {lh['1-60']['hazard']:.3f} over the first 60 days. By window (14 days / 28 days): {l14['1-7']['hazard']:.3f} / {l28['1-7']['hazard']:.3f} in days 1--7, {l14['8-14']['hazard']:.3f} / {l28['8-14']['hazard']:.3f} in days 8--14 and {l14['15-28']['hazard']:.3f} / {l28['15-28']['hazard']:.3f} in days 15--28. Difference between the windows (28 minus 14 days, bootstrap of papers): before the deadline, days 13 to 0 pooled, {100*hc['before_13_0']['h28']:.2f}\\% against {100*hc['before_13_0']['h14']:.2f}\\% per day (difference {100*hc['before_13_0']['diff']:.2f} points, 95\\% interval {100*hc['before_13_0']['lo']:.2f} to {100*hc['before_13_0']['hi']:.2f}); in the first week after it, {100*hc['after_1_7']['diff']:.1f} points ({100*hc['after_1_7']['lo']:.1f} to {100*hc['after_1_7']['hi']:.1f}).", size="\\scriptsize"))

T_rob = n
rob = R["robustness"]; AE = R["action_editors"]
def rr(lab, v):
    return f"{lab} & {num(v['W14']['n_papers'])} & {num(v['W28']['n_papers'])} & {ci(v['passthrough_q25'])} & {ci(v['passthrough_q50'])} & {ci(v['passthrough_q75'])} & {ci(v['passthrough_mean'])} & {pc(v['W14']['share_final48h'])} / {pc(v['W28']['share_final48h'])} & {pc(v['W14']['share_late'])} / {pc(v['W28']['share_late'])} \\\\"
base = dict(W14=G["W14"], W28=G["W28"], passthrough_q25=PT["q25"], passthrough_q50=PT["q50"], passthrough_q75=PT["q75"], passthrough_mean=PT["mean"])
body = [rr("Same-day panels (main sample)", base), rr("Last assignment within 2 days", rob["eps2"]), rr("Last assignment within 3 days", rob["eps3"]),
        rr("Under review in 2022--2023", rob["y2022_2023"]), rr("Under review in 2023--2025", rob["y2023_2025"]), rr("Under review in 2024--2025", rob["y2024_2025"]),
        rr("Panels with exactly three reviews", rob["exactly_three_reviews"]), rr("14-day panels labelled regular only", rob["W14_regular_label_only"])]
S.append(table("TMLR: pass-through of the window to the delivery date, alternative samples", "tab:tmlr-rob", "p{4.4cm}rrrrrrrr",
               "Sample & \\multicolumn{2}{c}{Papers} & \\multicolumn{4}{c}{Pass-through (95\\% interval)} & Final 48 h (\\%) & Late (\\%) \\\\\n & 14 d & 28 d & Lower quartile & Median & Upper quartile & Mean & 14 d / 28 d & 14 d / 28 d", body,
               "Pass-through: difference between the 28-day and the 14-day panels in the stated statistic of days since assignment, divided by 14; intervals from " + num(R["meta"]["nboot"]) + " bootstrap resamples of papers. In the second and third rows the deadline of a reviewer may precede the due date used by up to one and two days. Resampling action editors (" + str(AE['n_editors']) + ") instead of papers gives an interval of " + f"{AE['passthrough_q50_editor_bootstrap']['lo']:.2f}--{AE['passthrough_q50_editor_bootstrap']['hi']:.2f}" + " for the median of the main sample; within the " + str(AE['n_editors_both_windows']) + " action editors who handled panels of both kinds (" + str(AE['papers_in_both']['W14']) + " and " + str(AE['papers_in_both']['W28']) + " papers), the difference between their two medians, averaged with weights that reflect the number of panels, is " + f"{AE['within_editor_weighted_median_diff_over_14']:.2f}" + " of the 14 added days.",
               size="\\scriptsize", env="sidewaystable"))

T_len = n
fe = R["fe_log_chars"]; names = {"mid": "2--4 days before", "last48h": "24--48 hours before", "last24h": "Final 24 hours", "late": "Late"}
body = []
for b in ("mid", "last48h", "last24h", "late"):
    body.append(f"{names[b]} & " + " & ".join(f"{sgn(fe[k]['coef'][b]['b'],3)} ({fe[k]['coef'][b]['se']:.3f})" for k in ("all", "W14", "W28")) + " \\\\")
body.append("Reviews / papers & " + " & ".join(f"{num(fe[k]['n'])} / {num(fe[k]['papers'])}" for k in ("all", "W14", "W28")) + " \\\\")
body.append("Reviews by timing (early, 2--4 d, 24--48 h, final 24 h, late) & " + " & ".join(", ".join(str(fe[k]['bin_n'][b]) for b in ("early", "mid", "last48h", "last24h", "late")) for k in ("all", "W14", "W28")) + " \\\\")
lo = R["fe_log_chars_late_vs_on_time"]; f48 = R["fe_log_chars_final48h_vs_earlier_on_time"]; pw = R["fe_log_chars_per_week_late"]; ce = R["fe_claims_supported_late_vs_on_time"]; bw = R["log_chars_W28_minus_W14"]
body += ["\\midrule", "\\multicolumn{4}{@{}l}{\\emph{Single contrasts, both windows pooled}} \\\\",
         f"Late against on time & {sgn(lo['b'],3)} ({lo['se']:.3f}) & \\multicolumn{{2}}{{l}}{{{num(lo['n'])} reviews, {num(lo['papers'])} papers}} \\\\",
         f"Final 48 hours against earlier, on-time reviews only & {sgn(f48['b'],3)} ({f48['se']:.3f}) & \\multicolumn{{2}}{{l}}{{{num(f48['n'])} reviews, {num(f48['papers'])} papers}} \\\\",
         f"Per week of lateness & {sgn(pw['b'],3)} ({pw['se']:.3f}) & \\multicolumn{{2}}{{l}}{{{num(pw['n'])} reviews, {num(pw['papers'])} papers}} \\\\",
         f"Claims judged supported (0/1), late against on time & {sgn(ce['b'],3)} ({ce['se']:.3f}) & \\multicolumn{{2}}{{l}}{{{num(ce['n'])} reviews, {num(ce['papers'])} papers}} \\\\"]
S.append(table("TMLR: length of the report by timing, within manuscripts", "tab:tmlr-len", "p{6.2cm}rrr", "Timing of the review & Both windows & 14-day window & 28-day window", body,
               f"Dependent variable: logarithm of the number of characters of the review (all text fields). Manuscript fixed effects; reference category in the upper block: more than four days before the deadline; standard errors clustered by manuscript in parentheses. The last row uses the answer to the form question on whether the claims of the submission are supported by evidence ({pc(R['claims_supported_share'])}\\% answer yes). Between windows, reviews of the 28-day panels are {bw['diff']:.3f} log points longer than those of the 14-day panels (95\\% interval {bw['lo']:.3f} to {bw['hi']:.3f}); the manuscripts are longer too, so this difference is not an effect of the window.", size="\\scriptsize"))

T_cal = n
body = []
for r in C["scenarios"]:
    body.append(f"{r['eta']:.1f} & {r['L']:.1f} & {r['delta']} & {r['conference_profile']} & {r['journal_timing_conference_penalties']} & {r['journal_timing_journal_penalties']} & {r['journal_timing_no_penalty']} & {r['journal_timing_zero_hazard_beyond_27']} \\\\")
ch = C["check_profile_vs_observed"]
S.append(table("Journal scenarios with the delivery profile measured at TMLR instead of the conference profile: optimal window (days)", "tab:tmlr-cal", "rrrrrrrr",
               "$\\eta$ & $L$ & $\\delta$ & Conference profile & \\multicolumn{4}{c}{Delivery profile of TMLR} \\\\\n & & & (as in the manuscript) & conference shortfalls & journal shortfalls & no shortfall & no delivery $>$27 d ahead", body,
               f"Acceptance curve of the manuscript (two arms of the experiment of Chetty et al.\\ 2014) in every column. Delivery profile of TMLR: pooled daily hazards of the two windows in the last 14 days before the deadline, those of the 28-day panels for days 14--27, the mean hazard of days 14--27 ({pc(C['hazard_floor_beyond_27'],2)}\\% per day) beyond the data, and a late hazard of {C['mu_tmlr']:.3f} per day. Conference shortfalls: $\\Delta$ of the final 48 hours and of late reviews as in the manuscript; journal shortfalls: none in the final 48 hours and {C['Delta_late_tmlr']:.3f} for late reviews (Table~S{T_len}). The last column uses the journal shortfalls and sets the hazard to zero more than 27 days before the deadline. Windows searched from 3 to 60 days. The profile reproduces the two observed windows: late share {pc(ch['W14']['model_late'])}\\% (observed {pc(ch['W14']['observed_late'])}\\%) at 14 days and {pc(ch['W28']['model_late'])}\\% ({pc(ch['W28']['observed_late'])}\\%) at 28 days.", size="\\scriptsize"))

T_sel = n
SEL = R["selection"]; body = []
for c, lab in (("R", "Regular"), ("L", "Long")):
    for k, lab2 in (("same_day", "same-day panels"), ("other", "other panels")):
        d = SEL[c][k]
        body.append(f"{lab if k == 'same_day' else ''} & {lab2} & {num(d['n'])} & {pc(d['n']/SEL[c]['n_all_with_reviews'],1)} & {d['mean_reviews']:.2f} & {f1(d['first_review_median'])} & {f1(d['last_review_median'])} & {num(int(d['chars_median']))} & {pc(d['accept_share_of_decided'],1)} ({num(d['decided'])}) \\\\")
    ad = SEL[c]["accept_diff"]
    body.append(f" & \\multicolumn{{7}}{{l}}{{\\emph{{difference in acceptance, same-day minus other: {sgn(100*ad['diff'],1)} points (s.e.\\ {100*ad['se']:.1f})}}}} \\\\")
    if c == "R": body.append("\\midrule")
S.append(table("TMLR: same-day panels against the other papers with the same current length label", "tab:tmlr-sel", "llrrrrrrr",
               "Label & Panels & Papers & \\% of & Reviews & \\multicolumn{2}{c}{Median days to} & Median & Accepted, \\% \\\\\n & & & label & per paper & first review & last review & characters & (decided)", body,
               "Papers with a review due date, a length label and public reviews (" + num(A['with_due_and_reviews']) + "). Label: the current length label; regular submissions have a 14-day window and long ones a 28-day window. Same-day panels: the analysis sample with the stated current label (the " + str(G['W14_long_label']['n_papers']) + " 14-day panels relabelled long are in neither row of the lower block). Other panels: at least one reviewer assigned after the day the paper was put under review, or a due date set by hand. Days counted from the under-review date. Accepted: decided papers whose decision recommends acceptance, with or without revisions. Same-day panels are a larger share of regular than of long submissions, and among long submissions they are accepted less often than the others; the two arms of the comparison are therefore not selected in the same way.", size="\\scriptsize"))

sec = r"""
\section*{One journal, two windows: TMLR}
\begin{sloppypar}
\paragraph{Source.} Reviews: the snapshot of all public TMLR submissions and their replies taken on 13--14 October 2026 for a companion study, deposited at \url{https://github.com/rosadecsai/time-to-first-decision} (file \path{tfd_tmlr_all.jsonl.gz}; creation time of each review, signature, length in characters, short form answers; no text). Deadlines: the due date of the review invitation of each paper (\path{TMLR/Paper<n>/-/Review}), read from the OpenReview API on 27 November 2026 with the downloader of the main data (file \path{tmlr_review_duedates.csv}). These invitations are public for the papers put under review up to December 2025. In a sample of 1{,}223 reviews of 400 submissions the creation time used here is identical to the server-side timestamp \path{tcdate} in every case.

\paragraph{Rule.} In the workflow that OpenReview runs for the journal, the function that puts a paper under review posts the review invitation with a due date at 23:59:59 of the day two weeks later, or four weeks later when the submission is declared long, and the function that processes each reviewer assignment resets the due date in the same way and e-mails that date to the reviewer (\path{openreview-py}, module \path{journal}: \path{get_review_period_length}, \path{setup_under_review_submission} and \path{process/reviewer_assignment_process.py}; commit of 26 November 2026). The reviewer guide of the journal states the same rule. Table~S""" + str(T_acc) + r""" gives the construction of the sample that follows from it; three panels whose reviews are all dated 90 to 142 days before the recorded under-review date (the date is evidently a later one) are left out. The length label in the snapshot is the current one. In the 18 panels with a 14-day window and the label long whose edit history we read, the submission was declared regular when it was put under review and was relabelled long in a revision 24 to 87 days later; these panels are therefore regular submissions at the time of review, and they deliver like the others (Table~S""" + str(T_del) + r"""). We did not read the edit history of the remaining relabelled panels nor of the 28-day panels; a paper reviewed as regular, relabelled long later and whose last reviewer was assigned exactly 14 days after it was put under review would be counted as a 28-day panel.

\paragraph{Results.} Table~S""" + str(T_del) + r""" reports delivery by window, Table~S""" + str(T_haz) + r""" the daily hazards, Table~S""" + str(T_rob) + r""" the pass-through in alternative samples, Table~S""" + str(T_len) + r""" the length of the report by timing, Table~S""" + str(T_cal) + r""" the journal scenarios of the calibration recomputed with the delivery profile of TMLR, and Table~S""" + str(T_sel) + r""" how the same-day panels differ from the other papers. All are produced by \path{code/tmlr_analysis.py} and \path{code/tmlr_calibration.py}.
\end{sloppypar}

""" + "\n".join(S)
open(os.path.join(os.path.dirname(os.path.abspath(TAB)), "OnlineResource3_tmlr.tex"), "w").write(sec)
print("tables S%d-S%d" % (N0, n - 1), dict(acc=T_acc, deliv=T_del, haz=T_haz, rob=T_rob, length=T_len, cal=T_cal, sel=T_sel))
json.dump(dict(acc=T_acc, deliv=T_del, haz=T_haz, rob=T_rob, length=T_len, cal=T_cal, sel=T_sel), open(os.path.join(OUT, "tmlr_table_numbers.json"), "w"))
