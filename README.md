# The Optimal Review Deadline — code, data and results (ICLR 2018–2026)

Companion package for the manuscript on the optimal review window (theory + empirics
on time-stamped OpenReview reports). All numbers in the paper are produced by the
scripts in `code/` from the compact JSON files in `data/raw/` (one per ICLR year).

## Pipeline

```
python code/run_all.py 500        # data → E1–E7 → symbolic checks (≈ 12 min, 2 cores)
```

| Step | Script | Output |
| --- | --- | --- |
| Download (browser, logged-in OpenReview session) | `code/download_openreview_browser.js` | `data/raw/iclr_YYYY_notes.json[.gz]` |
| Build review-level dataset | `code/build_dataset.py` | `data/processed/reviews.csv`, `papers.csv`, `summary_by_year.csv` |
| E1 delivery distribution, E2 mixture, E3 Parkinson | `code/exp_timing.py` | Table 1, Table 2, Fig. 1, Fig. 2, `table_e3_*` |
| E4 quality by timing (manuscript FE), E5 windows, E6 tomography | `code/exp_quality.py` | Table 3, Table 4, Table 5, Fig. 3 |
| E7 calibration and T*, validation, Monte Carlo | `code/exp_model.py` | Table 6, Fig. 4, Fig. S1, `table_e7_*` |
| Propositions (SymPy + numerical) | `code/proofs_sympy.py` | 13/13 checks |
| Manuscript tables | `code/make_tables.py` | `manuscript/tables/*.tex` |

`code/config.py` holds the official windows (assignment → deadline, from iclr.cc reviewer
guides), the *effective* deadline dates detected from the hourly submission rate
(2022: 2 Nov 2021, one day after the published date; 2026: 31 Oct 2025 AoE, the
invitation duedate), the 48-hour atom definition and the seed (PCG64, 20260924).
`code/field_options.py` recovers long categorical answers (confidence, correctness)
that the compact downloader stored only as text metrics.

## Data

`data/raw/iclr_YYYY_notes.json.gz`: one record per note (reviews, meta-reviews,
decisions, comments) with `tcdate` (true creation time, not user-editable),
`cdate`, `tmdate`, `mdate`, anonymous signature, submission number, all short
fields (scores, confidence, flags) and, for long text fields,
`[chars, words, lines, bullets, question marks]`. No review text is stored.
Submissions carry their venue/venueid (decision, withdrawn, desk-rejected).

`data/processed/reviews.csv`: 209,389 official reviews (55,472 submissions):
timing relative to the effective deadline (`t_before`, days; `phase` = pre / atom
(last 48 h) / late; `bin` = early >4 d, mid 2–4 d, last48h, last24h, late), scores
(`score` = rating or recommendation, `confidence`, `soundness`, `presentation`,
`contribution`, `correctness`, novelty fields, ICLR-2020 self-reported effort),
text metrics (`words`, `bullets`, `qmarks`, `words_weak`, `words_quest`), panel
variables (`n_reviews`, `score_others`, `dev_panel`) and the decision.

## Main results (see `output/results.json`)

* 42–52 % of reviews arrive in the last 48 h before the deadline in every cohort;
  16–37 % arrive late. Quantiles of "days before the deadline" barely move with the
  window (median slope 0.07 d per day of window; 12–22 d): Parkinson's law.
* The daily hazard rises from below 1 %/day three weeks out to 59 % on the deadline day
  and is ≈0.29/day afterwards; a constant-hazard-from-assignment model is rejected
  in every cohort (log-likelihood ratio 0.28–0.67 per review).
* With manuscript fixed effects, last-24-hour reviews are 8.2 % shorter, list 0.54
  fewer bullet points, are 0.12 points less confident and 0.17 points more lenient,
  but do not disagree more with co-reviewers (κ by timing pair ≈ constant).
* Calibrated optimum: conferences (pre-committed reviewers) 14–25 days with
  compression (η = 1.3); journals (acceptance margin from Chetty, Saez & Sándor 2014)
  25–44 days; the loss from a too-short window is far larger than from a too-long one.

## Manuscript

`manuscript/deadline_paper.tex` (Springer sn-jnl, sn-basic author-year), `OnlineResource1.tex`, `OnlineResource2.tex`; compile with `pdflatex` + `bibtex` (class and .bst included).
