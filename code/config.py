"""Configuración del proyecto "plazo óptimo de revisión" (ICLR 2018-2026, OpenReview).

Ventanas oficiales de revisión (asignación -> plazo) tomadas de las guías de
revisores y CFP de iclr.cc (véase el plan de investigación). El plazo se fija a
las 23:59 AoE (UTC-12) del día oficial, es decir, 11:59 UTC del día siguiente.
"""
from __future__ import annotations

import datetime as dt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw"
OUT_DIR = ROOT / "data" / "processed"
TAB_DIR = ROOT / "output" / "tables"
FIG_DIR = ROOT / "output" / "figures"
for _d in (OUT_DIR, TAB_DIR, FIG_DIR):
    _d.mkdir(parents=True, exist_ok=True)

SEED = 20260924  # semilla PCG64 para bootstrap y Monte Carlo

# --- Ventanas oficiales -----------------------------------------------------
# assign: fecha de asignación de revisores (None si no consta en las páginas
# oficiales); deadline: día oficial del plazo (AoE); window: días.
YEARS = {
    2018: dict(assign=None, deadline=dt.date(2017, 11, 27), api="v1"),
    2019: dict(assign=None, deadline=dt.date(2018, 11, 2), api="v1"),
    2020: dict(assign=dt.date(2019, 10, 3), deadline=dt.date(2019, 10, 23), api="v1"),
    2021: dict(assign=dt.date(2020, 10, 12), deadline=dt.date(2020, 10, 28), api="v1"),
    2022: dict(assign=dt.date(2021, 10, 14), deadline=dt.date(2021, 11, 1), api="v1"),
    2023: dict(assign=dt.date(2022, 10, 12), deadline=dt.date(2022, 10, 24), api="v1"),
    2024: dict(assign=dt.date(2023, 10, 10), deadline=dt.date(2023, 10, 31), api="v2"),
    2025: dict(assign=dt.date(2024, 10, 14), deadline=dt.date(2024, 11, 3), api="v2"),
    2026: dict(assign=dt.date(2025, 10, 10), deadline=dt.date(2025, 11, 1), api="v2"),
}

AOE = dt.timezone(dt.timedelta(hours=-12))

# Plazo efectivo observado: día (AoE) cuyo final coincide con la caída brusca
# del ritmo horario de envíos. Coincide con el oficial salvo en:
#  - 2022: los informes se agrupan el 2-nov-2021, un día después del plazo de la
#    guía de revisores (extensión no documentada en la página oficial);
#  - 2026: la guía dice "November 01 2025", pero el duedate de la invitación
#    Official_Review es 2025-11-01T11:59Z = 31-oct-2025 23:59 AoE, y el ritmo
#    horario cae justo ahí.
# build_dataset.py comprueba que el ritmo horario cae en el instante fijado.
EFFECTIVE_DEADLINE = {y: v["deadline"] for y, v in YEARS.items()}
EFFECTIVE_DEADLINE[2022] = dt.date(2021, 11, 2)
EFFECTIVE_DEADLINE[2026] = dt.date(2025, 10, 31)


def _end_of_day_aoe(d: dt.date) -> dt.datetime:
    return dt.datetime(d.year, d.month, d.day, 23, 59, 59, tzinfo=AOE).astimezone(dt.timezone.utc)


def deadline_utc(year: int, effective: bool = True) -> dt.datetime:
    """Instante del plazo: 23:59:59 AoE del día del plazo (efectivo u oficial), en UTC."""
    d = EFFECTIVE_DEADLINE[year] if effective else YEARS[year]["deadline"]
    return _end_of_day_aoe(d)


# Átomo del plazo: las últimas ATOM_DAYS (48 h) antes de las 23:59 AoE. El día del
# plazo está en curso en alguna zona horaria durante 50 h (de las 00:00 en UTC+14
# a las 23:59 en UTC-12); las 48 h finales lo cubren salvo sus dos primeras horas.
ATOM_DAYS = 2.0


def deadline_date_window_utc(year: int) -> tuple[dt.datetime, dt.datetime]:
    """Intervalo UTC del átomo del plazo: [plazo - 48 h, plazo]."""
    end = _end_of_day_aoe(EFFECTIVE_DEADLINE[year])
    return end - dt.timedelta(days=ATOM_DAYS), end


def assign_utc(year: int) -> dt.datetime | None:
    """Inicio de la ventana: 00:00 AoE del día de asignación (None si no consta)."""
    a = YEARS[year]["assign"]
    if a is None:
        return None
    return dt.datetime(a.year, a.month, a.day, 0, 0, 0, tzinfo=AOE).astimezone(dt.timezone.utc)


def window_days(year: int, effective: bool = False) -> float | None:
    """Ventana (días) entre la asignación y el plazo publicado (o el efectivo)."""
    a = YEARS[year]["assign"]
    if a is None:
        return None
    d = EFFECTIVE_DEADLINE[year] if effective else YEARS[year]["deadline"]
    return float((d - a).days)


# Tramos temporales respecto al plazo (días antes del plazo, negativos = tarde)
BINS = [
    ("early", 4.0, float("inf")),       # > 4 días antes
    ("mid", 2.0, 4.0),                  # entre 2 y 4 días antes
    ("last48h", 1.0, 2.0),              # últimas 48 h (antes de las últimas 24)
    ("last24h", 0.0, 1.0),              # últimas 24 h
    ("late", float("-inf"), 0.0),       # tras el plazo
]
BIN_ORDER = [b[0] for b in BINS]

# Exclusión: revisores sustitutos ("substitute reviewing") de ICLR 2020, 23-oct a 4-nov
SUBSTITUTE_WINDOW = {2020: (dt.date(2019, 10, 24), dt.date(2019, 11, 5))}

# Ventanas de calibración externa (días) usadas en E7
CHETTY_DEADLINES = {"6 weeks": 45, "4 weeks": 28}
