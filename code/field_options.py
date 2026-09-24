"""Opciones de los campos categóricos largos (>120 caracteres) de los formularios
de revisión de ICLR, tomadas de las invitaciones Official_Review de la API
(consulta del 24-sep-2026). El descargador guarda para los textos largos solo
[caracteres, palabras, líneas, viñetas, interrogaciones]; con estas listas se
recupera el valor exacto por longitud y número de palabras.
"""

CONF_2018_2021 = [
    "5: The reviewer is absolutely certain that the evaluation is correct and very familiar with the relevant literature",
    "4: The reviewer is confident but not absolutely certain that the evaluation is correct",
    "3: The reviewer is fairly confident that the evaluation is correct",
    "2: The reviewer is willing to defend the evaluation, but it is quite likely that the reviewer did not understand central parts of the paper",
    "1: The reviewer's evaluation is an educated guess",
]

CONF_2022_2023 = [
    "1: You are unable to assess this paper and have alerted the ACs to seek an opinion from different reviewers.",
    "2: You are willing to defend your assessment, but it is quite likely that you did not understand the central parts of the submission or that you are unfamiliar with some pieces of related work. Math/other details were not carefully checked.",
    "3: You are fairly confident in your assessment. It is possible that you did not understand some parts of the submission or that you are unfamiliar with some pieces of related work. Math/other details were not carefully checked.",
    "4: You are confident in your assessment, but not absolutely certain. It is unlikely, but not impossible, that you did not understand some parts of the submission or that you are unfamiliar with some pieces of related work.",
    "5: You are absolutely certain about your assessment. You are very familiar with the related work and checked the math/other details carefully.",
]

CORRECTNESS_2022_2023 = [
    "1: The main claims of the paper are incorrect or not at all supported by theory or empirical results.",
    "2: Several of the paper’s claims are incorrect or not well-supported.",
    "3: Some of the paper’s claims have minor issues. A few statements are not well-supported, or require small changes to be made correct.",
    "4: All of the claims and statements are well-supported and correct.",
]

# año -> campo -> lista de opciones
OPTIONS = {
    2018: {"confidence": CONF_2018_2021},
    2019: {"confidence": CONF_2018_2021},
    2020: {},
    2021: {"confidence": CONF_2018_2021},
    2022: {"confidence": CONF_2022_2023, "correctness": CORRECTNESS_2022_2023},
    2023: {"confidence": CONF_2022_2023, "correctness": CORRECTNESS_2022_2023},
    2024: {"confidence": CONF_2022_2023},  # mismo texto en la API v2 (comprobado)
    2025: {},  # confidence numérica (1-5) en 2025-2026
    2026: {},
}


def _metrics(s: str):
    lines = [x.strip() for x in s.replace("\r\n", "\n").split("\n") if x.strip()]
    words = len(s.split())
    return (len(s), words)


def recover(year: int, field: str, metrics) -> str | None:
    """Devuelve la opción cuya (longitud, nº palabras) coincide con las métricas guardadas."""
    opts = OPTIONS.get(year, {}).get(field)
    if not opts or not metrics:
        return None
    key = (int(metrics[0]), int(metrics[1]))
    cands = [o for o in opts if _metrics(o) == key]
    if len(cands) == 1:
        return cands[0]
    cands = [o for o in opts if len(o) == key[0]]
    return cands[0] if len(cands) == 1 else None
