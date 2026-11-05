"""Set statistical symbols in figure text in italic (Nature style), added 3 Nov 2026 (Zira).

Imported by the *.py figure scripts. It patches matplotlib so that every text string
(titles, labels, legends, annotations) is passed through italicise() before drawing:
  p in 'Holm p', 'p = ...', 'p-values'  ->  italic capital P (Nature style: P values)
  n, D, V, f, z, q, T, k before =, <, >, <=, >=  ->  italic
  'ratio D', 'ln V', '(f;', 'V not tested', 'share f'  ->  italic
  Greek rho  ->  italic rho
Strings that already contain '$' (mathtext) are left alone. The font is fixed to DejaVu Sans,
with matching DejaVu italics for the symbols, so output matches across computers.
"""
import re

import matplotlib
import matplotlib.text as mtext

matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["mathtext.fontset"] = "dejavusans"
matplotlib.rcParams["mathtext.default"] = "it"

_RULES = [
    (re.compile(r"\bp-values?\b"), r"$P$ values"),
    (re.compile(r"\bHolm p\b"), r"Holm $P$"),
    (re.compile(r"(?<![\w$])p(?=\s*[=<>≤≥])"), r"$P$"),
    (re.compile(r"(?<![\w$])([nDVfzqTk])(?=\s*[=<>≤≥])"), r"$\1$"),
    (re.compile(r"\bratio D\b"), r"ratio $D$"),
    (re.compile(r"\bln V\b"), r"ln $V$"),
    (re.compile(r"\(f;"), r"($f$;"),
    (re.compile(r"\bshare f\b"), r"share $f$"),
    (re.compile(r"^V not tested$"), r"$V$ not tested"),
    (re.compile(r"ρ"), r"$\\rho$"),
    (re.compile(r"λ(\d)"), r"$\\lambda_\1$"),
]


def italicise(s):
    if not isinstance(s, str) or "$" in s or not s:
        return s
    for pat, rep in _RULES:
        s = pat.sub(rep, s)
    return s


_orig = mtext.Text.set_text


def _set_text(self, s):
    return _orig(self, italicise(s) if isinstance(s, str) else s)


mtext.Text.set_text = _set_text
