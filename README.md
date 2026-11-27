# The Optimal Review Deadline — code, data and results

Companion repository of "How Long Should Reviewers Get? Parkinson's Law and the Optimal
Review Deadline in 209,000 Timestamped Conference Reviews" (J. A. García, R. Rodriguez-Sánchez
and J. Fdez-Valdivia; *Scientometrics*, revised version R1, November 2026).

Two data sets, both read from the public OpenReview API and both without review text:

* **ICLR 2018–2026** (conference): 209,389 official reviews of 55,472 submissions
  (209,381 in the analysis sample).
* **TMLR 2022–2025** (journal): 17,280 reviews of 5,560 papers with a public review due date.
  The analysis uses the 1,115 "same-day panels" (3,426 reviews) in which all reviewers of a
  paper share one known window (14 days: 882 panels, 2,716 reviews; 28 days: 233 panels,
  710 reviews) and one known deadline.

Everything below runs with Python 3.10+ and `numpy`, `pandas`, `scipy`, `statsmodels`,
`matplotlib` and `sympy`. All random draws use fixed seeds (PCG64).

## 1. Original pipeline (ICLR)

```
python code/run_all.py 500        # data → E1–E7 → symbolic checks (≈ 12 min, 2 cores)
```

| Step | Script | Output |
| --- | --- | --- |
| Download (browser, signed-in OpenReview session) | `code/download_openreview_browser.js` | `data/raw/iclr_YYYY_notes.json.gz` |
| Build the review-level dataset | `code/build_dataset.py` | `data/processed/reviews.csv`, `papers.csv`, `summary_by_year.csv` |
| E1 delivery distribution, E2 mixture, E3 windows | `code/exp_timing.py` | delivery and mixture tables, Fig. 1, Fig. 2 |
| E4 visible effort by timing (manuscript fixed effects), E5 windows, E6 agreement | `code/exp_quality.py` | Tables 6–7 of the paper, Tables S1–S3 |
| E7 calibration, validation, Monte Carlo | `code/exp_model.py` | scenario table, Tables S4–S7, Fig. S1 |
| Propositions (SymPy and numerical) | `code/proofs_sympy.py` | 13/13 checks (four of the five comparative-statics signs are checked symbolically) |
| Manuscript tables | `code/make_tables.py` | LaTeX tables (the folder is given as an argument) |

`code/config.py` holds the official windows (assignment → deadline, from the iclr.cc reviewer
guides), the *effective* deadline dates detected from the hourly submission rate (2022: 2 Nov
2021, one day after the published date; 2026: 31 Oct 2025 AoE, the invitation due date), the
timing bins and the seed (20260924). `code/field_options.py` recovers long categorical answers
(confidence, correctness) that the compact downloader stored only as text metrics.

`output/results.json` and `output/figures/`, `output/tables/` are the archived outputs of this
pipeline. Figures 3 and 5 of the revised paper are redrawn from them by
`code/revised_figures.py` (section 2); `output/figures/fig4_optimal_deadline.*` is the figure of
the submitted version.

## 2. Analyses added in the revision

```
python code/revision_analyses.py . output/revision
python code/kappa_bootstrap.py . data/processed/reviews.csv.gz output/revision
python code/final_checks.py . output/revision data/processed/reviews.csv.gz
python code/revised_figures.py . output/revision
python code/tmlr_analysis.py data/tmlr/tfd_tmlr_all.jsonl.gz data/tmlr/tmlr_review_duedates.csv output/tmlr
python code/tmlr_calibration.py . output/tmlr
python code/tmlr_figure.py . output/tmlr
mkdir -p output/tex/tables
python code/make_tmlr_tables.py output/tmlr output/tex/tables 19
python code/make_rev_tables.py . output/revision output/tex/tables
```

About two minutes in all. Run in this order: `final_checks.py` and `revised_figures.py` use
the definitions of `revision_analyses.py`, and `make_rev_tables.py` reads the outputs of the
scripts above it.

| Script | What it does |
| --- | --- |
| `code/revision_analyses.py` | Cohort-level inference on the effect of the window (seven cohorts, t with 5 df, exact permutation, leave-one-out); exact reproduction of the 24 optima of the scenario table; alternative acceptance curves; partial pass-through and generalised closed form; alternative tails of the hazard; shortfall from 0 to 0.30; joint sensitivity (Latin hypercube, 20,000 draws, first-order indices, single-peak check); timing coefficients before and after 2024 and equivalence test; reviews per manuscript; Fig. 6. Computed from `output/results.json` and `data/processed/papers.csv.gz`. Seed 20261126. |
| `code/kappa_bootstrap.py` | Agreement (kappa) by timing pair with a bootstrap over manuscripts (4,000 resamples, seed 20261127), from the review-level file; it first reproduces the archived kappas exactly. |
| `code/final_checks.py` | Checks added after a last internal review: slope of the hazard across cohorts; slopes implied by full anchoring against observed slopes; pass-through for on-time reviews, all reviews and without the 2023 cohort; compression elasticity by contrast and the conference optimum as a function of it; published deadlines and windows instead of effective ones; fixed-effects regressions of scores by rating scale; shortfall capped at 1; functional form of the acceptance curve drawn in the joint analysis; late hazard of the reviewer model for β = 1, 0.9, 0.7. Tables S26–S29. |
| `code/revised_figures.py` | Figures 3 and 5 of the revised paper, redrawn from `output/results.json`. |
| `code/tmlr_analysis.py` | Same-day panels of TMLR: delivery by window; pass-through with a paper-clustered and an action-editor-clustered bootstrap (4,000 resamples) and within action editors; hazards by day and their contrasts; common-support test; report length with manuscript fixed effects; robustness; selection into the two windows. Seed 20261127. |
| `code/tmlr_calibration.py` | Journal scenarios with the delivery profile and the shortfall measured at TMLR. |
| `code/tmlr_figure.py` | Fig. 4. |
| `code/make_tmlr_tables.py`, `code/make_rev_tables.py` | Table 8 (TMLR), Table 5 (windows) and Table 10 (scenarios) of the paper, and Online Resource 3 (Tables S8–S29), written as LaTeX to `output/tex/`. The last argument of `make_tmlr_tables.py` is the number of the first TMLR table (S19). |
| `code/download_tmlr_duedates_browser.js` | The request that produced `data/tmlr/tmlr_review_duedates.csv` (ten requests; no reviewer identity is requested). |
| `code/audit_claims.py` | Recomputes the numbers of the revised text from the outputs and locates the sentence that states each one (365 numbers: 238 on the conference analyses, 127 on TMLR); checks that citations and bibliography match. It needs the LaTeX sources of the manuscript, which are not in this repository until publication: `python code/audit_claims.py . output/revision <manuscript folder> output/tmlr`. |

`output/revision/` and `output/tmlr/` hold the outputs of these scripts as used in the paper
(`*.json`, one CSV per analysis, the figures and the console logs). The per-review file of the
TMLR analysis sample is `output/tmlr/tmlr_same_day_panels_reviews.csv` (paper number, window,
days before the deadline, days since assignment, characters).

## 3. Data

`data/raw/iclr_YYYY_notes.json.gz`: one record per note (reviews, meta-reviews, decisions,
comments) with `tcdate` (creation time set by the server), `cdate`, `tmdate`, `mdate`, anonymous
signature, submission number, all short fields (scores, confidence, flags) and, for long text
fields, `[chars, words, lines, bullets, question marks]`. No review text is stored.

`data/processed/reviews.csv.gz`: 209,389 official reviews: timing relative to the effective
deadline (`t_before`, days) and to the published one (`t_before_official`); `bin` = early (four
days or more before the deadline), mid (2–4 days), last48h (24–48 hours before), last24h, late;
`phase` = pre / atom (final 48 hours) / late; scores, text metrics (`words`, `bullets`,
`qmarks`, `words_weak`, `words_quest`), panel variables (`n_reviews`, `score_others`, `dev_panel`), `in_sample` and the
decision. `papers.csv.gz`, `summary_by_year.csv`: paper- and cohort-level files.

`data/tmlr/tfd_tmlr_all.jsonl.gz`: snapshot of all public TMLR submissions and their replies
(13–14 October 2026), copied from `rosadecsai/time-to-first-decision`. One line per submission;
reviews carry their creation time, the anonymous per-paper signature, length and short form
answers. No text.

`data/tmlr/tmlr_review_duedates.csv`: due date of the review invitation of each paper
(`TMLR/Paper<n>/-/Review`), read on 27 November 2026: 5,651 papers put under review up to
December 2025. Columns: `number`, `due_ms` (UTC, milliseconds), `end_of_day` (1 when the date is
23:59:59 UTC, the value the journal workflow sets; 110 papers with reviews carry a date set by
hand and are left out).

How the TMLR due date is set, and why the same-day panels identify the window: the OpenReview
journal workflow (`openreview-py`, module `journal`) posts the review invitation with a due date
at 23:59:59 of the day two weeks after the paper is put under review, four weeks for a long
submission, and resets it in the same way at every reviewer assignment. The invitation keeps
the date of the last assignment, so a paper whose due date falls 14 (28) days after the day it
was put under review had all its reviewers assigned on that day. Three panels with a review
dated before the paper was put under review are dropped.

## 4. Main results

* **ICLR, timing.** 42–52 % of reviews arrive in the final 48 hours before the deadline in
  every cohort and 16–37 % arrive late. The daily hazard of delivery rises from below 1 % two
  to three weeks before the deadline to 59 % in the final 24 hours. Across seven cohorts with
  windows of 12–21 days (not randomised) the median delivery date of on-time reviews moves
  0.93 days per day of window (cohort-level 95 % interval 0.85–1.00; 0.94 with late reviews
  included; 0.89, 0.77–1.01, with the published deadlines).
* **TMLR, timing.** Reviewers given 28 days instead of 14 deliver 13.8 days later at the median
  (0.99 days per day of window; 95 % interval 0.95–1.04, the same when resampling action
  editors). Before the deadline the daily hazards of the two windows coincide when aligned on
  it (3.4 % against 3.2 % per day); after it they differ (16.5 % against 10.6 % on the first
  late day). The deadline is softer than at the conference: 15–19 % of reviews arrive in the
  final 48 hours and 53–57 % after the deadline. The windows are assigned by manuscript
  length, not at random, so the slope bounds the effect of the window from above.
* **Visible effort.** Within ICLR manuscripts, final-24-hour reviews are 8.2 % shorter, less
  confident and more lenient (0.13 points on the 1–10 scale); their scores deviate from their
  co-reviewers' no more than early ones, and pairs of late reviews agree slightly less than
  pairs of early ones (kappa lower by 0.05; 0.02–0.09). At TMLR last-minute reviews are not
  shorter; late ones are.
* **Calibration (conditional on stated inputs).** For conferences the window the model favours
  depends on how fast the loss from last-minute reviews grows as windows shorten, an added
  assumption that the data inform weakly: at the central delay cost, 5 days if it does not
  grow, 14 if it grows as the comparison with the 2020–2022 cohorts suggests and 19 if it grows
  as the comparison with the 2024–2026 cohorts suggests. For journals that recruit reviewers
  the result is a rule, not a number: the illustrative range (23–44 days; 25–49 with the
  delivery profile of TMLR) rests on an acceptance curve known at two windows from one
  experiment (Chetty, Saez and Sándor 2014), and with the functional form of that curve drawn
  along with the other inputs the central 90 % of the optima runs from 21 to 90 days.

Code: MIT licence. Data derived from the public OpenReview API.
