"""Construye el conjunto de datos agregado por informe a partir de los JSON
compactos descargados de OpenReview (iclr_YYYY_notes.json[.gz]).

Salida (data/processed/):
  reviews.csv  una fila por Official_Review (sin texto)
  papers.csv   una fila por envío (decisión, nº de informes, retirada)
  summary_by_year.csv  recuento y comprobaciones por año
"""
from __future__ import annotations

import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import field_options as FO  # noqa: E402

NUM_FIELDS = [
    "rating", "recommendation", "confidence", "correctness",
    "technical_novelty_and_significance", "empirical_novelty_and_significance",
    "soundness", "presentation", "contribution",
]

# Campos de texto largo por año (para longitud total del informe)
TEXT_FIELDS = {
    "review", "main_review", "summary_of_the_paper", "summary_of_the_review",
    "strength_and_weaknesses", "clarity,_quality,_novelty_and_reproducibility",
    "summary", "strengths", "weaknesses", "questions", "details_of_ethics_concerns",
    "limitations", "review_summary",
}
WEAK_FIELDS = {"weaknesses", "strength_and_weaknesses"}
QUEST_FIELDS = {"questions"}

# ICLR 2020: escalas ordinales de esfuerzo autodeclarado
EXP_2020 = {
    "I do not know much about this area.": 1,
    "I have read many papers in this area.": 2,
    "I have published one or two papers in this area.": 3,
    "I have published in this field for several years.": 4,
}
THOR_2020 = {
    "I made a quick assessment of this paper.": 1,
    "I read the paper thoroughly.": 2,
    "I read the paper at least twice and used my best judgement in assessing the paper.": 3,
}
DERIV_2020 = {
    "N/A": np.nan,
    "I did not assess the derivations or theory.": 0,
    "I assessed the sensibility of the derivations and theory.": 1,
    "I carefully checked the derivations and theory.": 2,
}
EXPER_2020 = {
    "N/A": np.nan,
    "I did not assess the experiments.": 0,
    "I assessed the sensibility of the experiments.": 1,
    "I carefully checked the experiments.": 2,
}


def load_year(year: int) -> dict:
    for name in (f"iclr_{year}_notes.json.gz", f"iclr_{year}_notes.json"):
        p = C.RAW_DIR / name
        if p.exists():
            if name.endswith(".gz"):
                with gzip.open(p, "rt", encoding="utf-8") as f:
                    return json.load(f)
            return json.loads(p.read_text(encoding="utf-8"))
    raise FileNotFoundError(f"no hay datos para {year}")


def num_prefix(v) -> float:
    """'6: marginally above...' -> 6 ; '3 good' -> 3 ; 5 -> 5."""
    if v is None:
        return np.nan
    if isinstance(v, (int, float)):
        return float(v)
    m = re.match(r"\s*(-?\d+(?:\.\d+)?)", str(v))
    return float(m.group(1)) if m else np.nan


def decision_code(s) -> str | None:
    """Normaliza la decisión a accept / reject / withdrawn / desk / other."""
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return None
    t = str(s).lower()
    if "withdraw" in t:
        return "withdrawn"
    if "desk" in t:
        return "desk_rejected"
    if "accept" in t or "oral" in t or "poster" in t or "spotlight" in t or "invite to workshop" in t:
        return "accept" if "workshop" not in t else "workshop"
    if "reject" in t:
        return "reject"
    if "submitted to" in t:
        return "pending"
    return "other"


def build_year(year: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    d = load_year(year)
    notes = d["notes"]
    subs = d["subs"]
    dl = C.deadline_utc(year)                       # plazo efectivo (23:59 AoE)
    dl_off = C.deadline_utc(year, effective=False)  # plazo oficial publicado
    day_start, day_end = C.deadline_date_window_utc(year)
    asg = C.assign_utc(year)
    dl_ms = dl.timestamp() * 1000
    dl_off_ms = dl_off.timestamp() * 1000
    day_start_ms = day_start.timestamp() * 1000
    asg_ms = asg.timestamp() * 1000 if asg else np.nan

    # --- envíos ---
    prow = {}
    for s in subs:
        vs = s.get("vs") or []
        prow[s["i"]] = dict(
            year=year, forum=s["i"], paper_num=s.get("n"), sub_cdate=s.get("c"), sub_tcdate=s.get("tc"),
            withdrawn=any("Withdrawn" in v for v in vs) or str(s.get("venueid", "")).endswith("Withdrawn_Submission"),
            desk_rejected=any("Desk_Rejected" in v for v in vs) or str(s.get("venueid", "")).endswith("Desk_Rejected_Submission"),
            venue=s.get("venue"), venueid=s.get("venueid"), area=s.get("area"), n_direct_replies=s.get("nrep"),
        )
    # decisiones y meta-reviews
    dec, meta_tc, withdraw_tc = {}, {}, {}
    for n in notes:
        v = n["v"]
        if v == "Decision":
            dec[n["f"]] = (n.get("S") or {}).get("decision")
        elif v == "Meta_Review" and n["f"] not in dec:
            S = n.get("S") or {}
            if "recommendation" in S:
                dec[n["f"]] = S["recommendation"]
            meta_tc[n["f"]] = n["tc"]
        elif v in ("Withdraw", "Withdraw_Submission", "Withdrawal"):
            withdraw_tc[n["f"]] = min(n["tc"], withdraw_tc.get(n["f"], n["tc"]))
    for f, p in prow.items():
        raw = dec.get(f)
        if raw is None and p["venue"]:
            raw = p["venue"]
        p["decision_raw"] = raw
        code = decision_code(raw)
        if p["withdrawn"]:
            code = "withdrawn"
        elif p["desk_rejected"]:
            code = "desk_rejected"
        p["decision"] = code
        p["withdraw_tc"] = withdraw_tc.get(f)
        p["withdrawn_before_assign"] = bool(withdraw_tc.get(f) and asg and withdraw_tc[f] < asg_ms)
        p["withdrawn_before_deadline"] = bool(withdraw_tc.get(f) and withdraw_tc[f] < dl_ms)

    # --- informes ---
    rows = []
    for n in notes:
        if n["v"] != "Official_Review":
            continue
        S = n.get("S") or {}
        L = n.get("L") or {}
        r = dict(
            year=year, review_id=n["i"], forum=n["f"], paper_num=n.get("n"), reviewer=n.get("s"),
            tcdate=n["tc"], cdate=n.get("c"), tmdate=n.get("tm"), mdate=n.get("md"),
        )
        r["t_before"] = (dl_ms - n["tc"]) / 86400e3            # días antes del plazo efectivo (neg. = tarde)
        r["t_before_official"] = (dl_off_ms - n["tc"]) / 86400e3
        r["t_since_assign"] = (n["tc"] - asg_ms) / 86400e3 if asg else np.nan
        r["late"] = r["t_before"] < 0
        # átomo del plazo: últimas 48 h antes de las 23:59 AoE (cubre el día del plazo en casi toda zona horaria)
        r["on_deadline_date"] = (n["tc"] >= day_start_ms) and (n["tc"] <= dl_ms)
        r["before_deadline_date"] = n["tc"] < day_start_ms
        r["edited"] = (n.get("tm") or n["tc"]) - n["tc"] > 3600e3
        r["edit_delay_days"] = ((n.get("tm") or n["tc"]) - n["tc"]) / 86400e3
        # campos numéricos (prefijo numérico); recuperar los largos desde L
        for k in NUM_FIELDS:
            val = S.get(k)
            if val is None and k in L:
                val = FO.recover(year, k, L[k])
            r[k] = num_prefix(val) if val is not None else np.nan
        # esfuerzo autodeclarado ICLR 2020
        if year == 2020:
            r["experience"] = EXP_2020.get(S.get("experience_assessment"), np.nan)
            r["thoroughness"] = THOR_2020.get(S.get("review_assessment:_thoroughness_in_paper_reading"), np.nan)
            r["check_deriv"] = DERIV_2020.get(S.get("review_assessment:_checking_correctness_of_derivations_and_theory"), np.nan)
            r["check_exper"] = EXPER_2020.get(S.get("review_assessment:_checking_correctness_of_experiments"), np.nan)
        # métricas de texto
        chars = words = lines = bullets = qmarks = 0
        w_weak = w_quest = 0
        for k in TEXT_FIELDS:
            if k in L:
                c_, w_, l_, b_, q_ = L[k]
            elif k in S and isinstance(S[k], str):
                s_ = S[k]
                c_, w_, l_, b_, q_ = len(s_), len(s_.split()), 1, 0, s_.count("?")
            else:
                continue
            chars += c_; words += w_; lines += l_; bullets += b_; qmarks += q_
            if k in WEAK_FIELDS:
                w_weak += w_
            if k in QUEST_FIELDS:
                w_quest += w_
        r.update(chars=chars, words=words, lines=lines, bullets=bullets, qmarks=qmarks,
                 words_weak=w_weak, words_quest=w_quest)
        eth = S.get("flag_for_ethics_review")
        r["ethics_flag"] = (not str(eth).upper().startswith("['NO") and not str(eth).upper().startswith("NO")) if eth is not None else np.nan
        rows.append(r)
    rev = pd.DataFrame(rows)
    pap = pd.DataFrame(list(prow.values()))

    # --- variables de panel ---
    rev["score"] = rev["rating"].where(rev["rating"].notna(), rev["recommendation"])
    g = rev.groupby("forum")["score"]
    rev["n_reviews"] = g.transform("count")
    rev["score_others"] = (g.transform("sum") - rev["score"]) / (rev["n_reviews"] - 1)
    rev["dev_panel"] = (rev["score"] - rev["score_others"]).abs()
    rev = rev.merge(pap[["forum", "decision", "withdrawn", "desk_rejected", "withdrawn_before_assign", "withdrawn_before_deadline", "withdraw_tc"]], on="forum", how="left")
    # tramo temporal
    bins = pd.Series(pd.NA, index=rev.index, dtype="object")
    for name, lo, hi in C.BINS:
        bins[(rev["t_before"] >= lo) & (rev["t_before"] < hi)] = name
    rev["bin"] = pd.Categorical(bins, categories=C.BIN_ORDER)
    # exclusiones
    rev["substitute"] = False
    if year in C.SUBSTITUTE_WINDOW:
        a, b = C.SUBSTITUTE_WINDOW[year]
        t = pd.to_datetime(rev["tcdate"], unit="ms", utc=True).dt.date
        rev["substitute"] = (t >= a) & (t <= b)
    rev["in_sample"] = ~rev["withdrawn_before_assign"].fillna(False).astype(bool)
    # fase respecto al plazo efectivo: pre (antes de las últimas 48 h), atom (últimas 48 h), late (tras las 23:59 AoE)
    rev["phase"] = np.where(rev["late"], "late", np.where(rev["on_deadline_date"], "atom", "pre"))
    pap["n_reviews"] = pap["forum"].map(rev.groupby("forum").size()).fillna(0).astype(int)
    pap["mean_score"] = pap["forum"].map(rev.groupby("forum")["score"].mean())
    return rev, pap


def main(years=None):
    years = years or [y for y in C.YEARS if any((C.RAW_DIR / f"iclr_{y}_notes.json{ext}").exists() for ext in (".gz", ""))]
    R, P, summ = [], [], []
    for y in years:
        rev, pap = build_year(y)
        R.append(rev); P.append(pap)
        dl = C.deadline_utc(y)
        # comprobaciones: (i) el ritmo horario en el átomo (48 h) supera con claridad
        # al de las 24 h siguientes; (ii) ritmo en
        # las 6 h anteriores al instante AoE frente a las 6 h posteriores (informativo:
        # en 2020 no hay caída en el instante AoE, solo en el día).
        hrs = np.floor(-rev["t_before"] * 24).astype(int)   # horas desde el plazo (neg. = antes)
        rate_day = float(rev["on_deadline_date"].sum()) / (24.0 * C.ATOM_DAYS)
        rate_after = float(((hrs >= 0) & (hrs < 24)).sum()) / 24.0
        cliff_day = rate_day / max(rate_after, 1e-9)
        before6 = int(((hrs >= -6) & (hrs < 0)).sum()) / 6.0
        after6 = int(((hrs >= 0) & (hrs < 6)).sum()) / 6.0
        cliff = before6 / max(after6, 1e-9)
        assert cliff_day > 1.3, (y, "el día del plazo efectivo no destaca", rate_day, rate_after)
        days = pd.to_datetime(rev["tcdate"], unit="ms", utc=True).dt.floor("D")
        peak_day = days.value_counts().idxmax().date()
        summ.append(dict(
            year=y, n_papers=len(pap), n_reviews=len(rev), n_in_sample=int(rev["in_sample"].sum()),
            window_days=C.window_days(y), window_days_effective=C.window_days(y, effective=True),
            official_deadline=str(C.YEARS[y]["deadline"]),
            effective_deadline=str(C.EFFECTIVE_DEADLINE[y]), peak_day_utc=str(peak_day),
            cliff_day_ratio=round(cliff_day, 2), cliff_hour_ratio=round(cliff, 2),
            share_deadline_date=float(rev["on_deadline_date"].mean()),
            share_last24h=float((rev["bin"] == "last24h").mean()),
            share_last48h=float(rev["bin"].isin(["last24h", "last48h"]).mean()),
            share_late=float(rev["late"].mean()),
            median_t_before=float(rev["t_before"].median()),
            share_score_known=float(rev["score"].notna().mean()), share_conf_known=float(rev["confidence"].notna().mean()),
            mean_words=float(rev["words"].mean()),
        ))
        s = summ[-1]
        print(f"{y}: {len(pap)} envíos, {len(rev)} informes, pico {peak_day} (plazo efectivo {C.EFFECTIVE_DEADLINE[y]}), "
              f"día del plazo {s['share_deadline_date']:.1%}, últimas 48 h {s['share_last48h']:.1%}, tardíos {s['share_late']:.1%}, "
              f"conf. conocida {s['share_conf_known']:.1%}")
    rev = pd.concat(R, ignore_index=True)
    pap = pd.concat(P, ignore_index=True)
    rev.to_csv(C.OUT_DIR / "reviews.csv", index=False)
    pap.to_csv(C.OUT_DIR / "papers.csv", index=False)
    pd.DataFrame(summ).to_csv(C.OUT_DIR / "summary_by_year.csv", index=False)
    return rev, pap


if __name__ == "__main__":
    main()
