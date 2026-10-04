"""Download GEO series matrices (GSE20916, GSE37364) and the GPL570 annotation into data/raw/."""
import os

import GEOparse
import requests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")

for gse in ("GSE20916", "GSE37364"):
    dst = os.path.join(RAW, f"{gse}_series_matrix.txt.gz")
    if os.path.exists(dst):
        print(f"{gse}: already present ({os.path.getsize(dst) / 1e6:.1f} MB)")
        continue
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}/matrix/"
           f"{gse}_series_matrix.txt.gz")
    r = requests.get(url, timeout=600)
    r.raise_for_status()
    with open(dst, "wb") as fh:
        fh.write(r.content)
    print(f"{gse}: downloaded {len(r.content) / 1e6:.1f} MB from {url}")

gpl = GEOparse.get_GEO(geo="GPL570", destdir=RAW, silent=True)
print(f"GPL570: {len(gpl.table)} probes; columns include {list(gpl.table.columns)[:12]}")
for f in sorted(os.listdir(RAW)):
    print(f"  {f:40s} {os.path.getsize(os.path.join(RAW, f)) / 1e6:7.1f} MB")
