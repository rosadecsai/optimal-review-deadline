"""Utilidades comunes: carga de datos, bootstrap por manuscrito, regresiones con
efectos fijos de manuscrito y errores agrupados, salida de tablas."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402

RESULTS = C.ROOT / "output" / "results.json"


def load_reviews(in_sample: bool = True) -> pd.DataFrame:
    df = pd.read_csv(C.OUT_DIR / "reviews.csv", low_memory=False)
    df["bin"] = pd.Categorical(df["bin"], categories=C.BIN_ORDER)
    df["phase"] = pd.Categorical(df["phase"], categories=["pre", "atom", "late"])
    if in_sample:
        df = df[df["in_sample"].astype(bool)].copy()
    df["log_words"] = np.log(df["words"].clip(lower=1))
    return df


def load_papers() -> pd.DataFrame:
    return pd.read_csv(C.OUT_DIR / "papers.csv", low_memory=False)


def rng(seed: int | None = None) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(C.SEED if seed is None else seed))


def cluster_bootstrap(df: pd.DataFrame, stat, n_boot: int = 500, cluster: str = "forum", seed: int | None = None):
    """Bootstrap por conglomerados (manuscritos). stat(df) -> escalar o array.
    Devuelve (estimación, IC 2.5 %, IC 97.5 %)."""
    g = rng(seed)
    ids = df[cluster].unique()
    groups = {k: np.asarray(v) for k, v in df.groupby(cluster).groups.items()}
    est = np.asarray(stat(df), dtype=float)
    reps = []
    for _ in range(n_boot):
        pick = g.choice(ids, size=len(ids), replace=True)
        idx = np.concatenate([groups[k] for k in pick])
        reps.append(np.asarray(stat(df.loc[idx]), dtype=float))
    reps = np.vstack(reps)
    lo, hi = np.nanpercentile(reps, [2.5, 97.5], axis=0)
    return est, lo, hi


def fe_regression(df: pd.DataFrame, y: str, xcols: list[str], fe: str = "forum", cluster: str | None = "forum",
                  weights: str | None = None):
    """OLS de y sobre xcols con efectos fijos de `fe` (transformación within) y
    errores estándar agrupados por `cluster`. Devuelve DataFrame coef/se/t/p y n."""
    import statsmodels.api as sm

    d = df[[y, fe] + xcols].dropna().copy()
    # solo grupos con >= 2 observaciones aportan identificación within
    cnt = d.groupby(fe)[y].transform("size")
    d = d[cnt >= 2]
    if d.empty:
        return None
    Y = d[y] - d.groupby(fe)[y].transform("mean")
    X = d[xcols] - d.groupby(fe)[xcols].transform("mean")
    X = X.loc[:, (X.abs() > 1e-12).any(axis=0)]  # quitar columnas sin variación within
    model = sm.OLS(Y.values, X.values)
    if cluster:
        res = model.fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(d[cluster])[0]})
    else:
        res = model.fit(cov_type="HC1")
    # corrección de grados de libertad por los efectos fijos absorbidos
    n_groups = d[fe].nunique()
    k = X.shape[1]
    dof_adj = (len(d) - k) / (len(d) - k - n_groups + 1)
    se = res.bse * np.sqrt(dof_adj)
    out = pd.DataFrame({"coef": res.params, "se": se, "t": res.params / se}, index=X.columns)
    from scipy import stats as st
    out["p"] = 2 * st.norm.sf(out["t"].abs())
    out.attrs["n"] = len(d)
    out.attrs["n_groups"] = n_groups
    out.attrs["r2_within"] = res.rsquared
    return out


def save_result(key: str, value):
    """Guarda resultados clave en output/results.json (para el manuscrito)."""
    data = {}
    if RESULTS.exists():
        data = json.loads(RESULTS.read_text())
    data[key] = value
    RESULTS.write_text(json.dumps(data, indent=1, default=_json_default))


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, pd.DataFrame):
        return o.to_dict(orient="records")
    return str(o)


def to_latex(df: pd.DataFrame, path: Path, caption: str, label: str, float_fmt: str = "%.2f", index: bool = False):
    """Tabla LaTeX (booktabs) junto con CSV."""
    df.to_csv(path.with_suffix(".csv"), index=index)
    tex = df.to_latex(index=index, float_format=float_fmt, escape=True, caption=caption, label=label)
    tex = tex.replace("\\toprule", "\\toprule").replace("\\hline", "")
    path.with_suffix(".tex").write_text(tex)
