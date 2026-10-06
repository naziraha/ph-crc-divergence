"""Fetch GEO series + sample METADATA only (brief SOFT view, no expression values).

Usage: python scripts/fetch_geo_metadata.py GSE117606 GSE41657
Writes data/raw/<GSE>_series_brief.txt and data/raw/<GSE>_samples_brief.txt and prints a summary.
"""
import os
import re
import sys
from collections import Counter

import requests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"


def get(acc, targ):
    r = requests.get(URL, params={"acc": acc, "targ": targ, "form": "text", "view": "brief"}, timeout=120)
    r.raise_for_status()
    return r.text


def parse_samples(text):
    samples, cur = [], None
    for line in text.splitlines():
        if line.startswith("^SAMPLE"):
            cur = {"gsm": line.split("=")[1].strip(), "ch": []}
            samples.append(cur)
        elif cur is not None and line.startswith("!Sample_"):
            k, _, v = line.partition(" = ")
            k = k.strip()
            if k == "!Sample_characteristics_ch1":
                cur["ch"].append(v.strip())
            elif k in ("!Sample_title", "!Sample_source_name_ch1", "!Sample_platform_id", "!Sample_description",
                       "!Sample_relation"):
                cur.setdefault(k, []).append(v.strip())
    return samples


for gse in sys.argv[1:]:
    s_txt, g_txt = get(gse, "self"), get(gse, "gsm")
    open(os.path.join(RAW, f"{gse}_series_brief.txt"), "w", encoding="utf-8").write(s_txt)
    open(os.path.join(RAW, f"{gse}_samples_brief.txt"), "w", encoding="utf-8").write(g_txt)
    print("=" * 80)
    for key in ("!Series_title", "!Series_platform_id", "!Series_summary", "!Series_overall_design",
                "!Series_relation", "!Series_contributor", "!Series_pubmed_id", "!Series_submission_date"):
        vals = [l.partition(" = ")[2] for l in s_txt.splitlines() if l.startswith(key)]
        if vals:
            print(f"{key}: " + " | ".join(v[:600] for v in vals))
    sm = parse_samples(g_txt)
    print(f"samples: {len(sm)}")
    print("platform:", Counter(x.get("!Sample_platform_id", ["?"])[0] for x in sm))
    print("source:", Counter(x.get("!Sample_source_name_ch1", ["?"])[0] for x in sm).most_common(10))
    keys = Counter()
    for x in sm:
        for c in x["ch"]:
            keys[c.split(":")[0].strip().lower()] += 1
    print("characteristic fields:", dict(keys))
    for fld in keys:
        vals = Counter(next((c.split(":", 1)[1].strip() for c in x["ch"] if c.split(":")[0].strip().lower() == fld), "")
                       for x in sm)
        if len(vals) <= 20:
            print(f"  {fld}: {dict(vals)}")
        else:
            print(f"  {fld}: {len(vals)} distinct, e.g. {list(vals)[:8]}")
    print("titles (first 12):", [x.get("!Sample_title", [""])[0] for x in sm][:12])
    rel = Counter(r for x in sm for r in x.get("!Sample_relation", []) if not r.startswith(("SRA", "BioSample")))
    if rel:
        print("sample relations:", dict(rel.most_common(8)))
