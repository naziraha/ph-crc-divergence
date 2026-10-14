"""
Download every public input file into data/raw/ and compare it with the SHA-256 recorded in data/PROVENANCE.md.

Sources
  GEO series matrices   GSE41258, GSE20916, GSE37364, GSE117606 (NCBI GEO FTP)
  GEO platform tables   GPL96, GPL570, GPL25373 (GEOparse, as in the analysis scripts)
  NCBI gene_info        Homo_sapiens.gene_info.gz (Entrez -> symbol, HGNC nomenclature-authority symbols)
  UCSC Xena TCGA hub    TCGA-COAD HiSeqV2, clinical matrix, MC3 gene-level mutations, GISTIC2 gene-level copy number
  data/raw/GSE117606_cel_headers.csv is shipped with the repository (derived from the GEO CEL archive by
  scripts/gse117606_cel_headers.py) and is only checked here.

Checksums are those of the files downloaded for the paper (dates in data/PROVENANCE.md).
  - NCBI gene_info (updated daily) and the Xena clinical matrix (updated occasionally) are expected to change:
    a difference is reported as a WARNING with an explanation, and the run continues.
  - Any other difference is an ERROR (exit code 1), because the input is then not the one used for the paper;
    --allow-mismatch continues anyway.
The exact files used for the paper are archived with the code release (see README, "Data").

Usage: python data/download_data.py [--allow-mismatch] [--check-only]
"""
import argparse
import gzip
import hashlib
import os
import re
import shutil
import sys

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")
GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}nnn/{gse}/matrix/{gse}_series_matrix.txt.gz"
XENA = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
FILES = {  # local name -> (url, decompress gzip after download)
    **{f"{g}_series_matrix.txt.gz": (GEO_FTP.format(stub=g[:-3], gse=g), False)
       for g in ("GSE41258", "GSE20916", "GSE37364", "GSE117606")},
    "Homo_sapiens.gene_info.gz": ("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz", False),
    "TCGA/TCGA.COAD.sampleMap_HiSeqV2.gz": (XENA + "TCGA.COAD.sampleMap%2FHiSeqV2.gz", False),
    "TCGA/TCGA.COAD.sampleMap_COAD_clinicalMatrix.txt": (XENA + "TCGA.COAD.sampleMap%2FCOAD_clinicalMatrix", False),
    "TCGA/COAD_mc3_gene_level.txt": (XENA + "mc3_gene_level%2FCOAD_mc3_gene_level.txt.gz", True),
    "TCGA/TCGA.COAD.sampleMap_Gistic2_all_data_by_genes.gz": (XENA + "TCGA.COAD.sampleMap%2FGistic2_CopyNumber_Gistic2_all_data_by_genes.gz", False),
}
PLATFORMS = ("GPL96", "GPL570", "GPL25373")
VOLATILE = {  # sources that are updated in place; a newer version is expected and only warned about
    "GPL96.txt": "GEO appends every newly submitted sample and series that uses the platform to the header of the "
                "platform record, so the file changes over time. The probe annotation table (ID, Gene Symbol) is "
                "what the pipeline reads; for GPL570 it was checked to be identical to the version of 14 October 2026.",
    "GPL570.txt": "GEO appends every newly submitted sample and series that uses the platform to the header of the "
                "platform record, so the file changes over time. The probe annotation table (ID, Gene Symbol) is "
                "what the pipeline reads; for GPL570 it was checked to be identical to the version of 14 October 2026.",
    "GPL25373.txt": "GEO appends every newly submitted sample and series that uses the platform to the header of the "
                "platform record, so the file changes over time. The probe annotation table (ID, Gene Symbol) is "
                "what the pipeline reads; for GPL570 it was checked to be identical to the version of 14 October 2026.",
    "Homo_sapiens.gene_info.gz": "NCBI regenerates gene_info daily (new genes, renamed symbols, retired IDs). A newer "
                                 "file can change the Entrez-to-symbol mapping of GSE117606 and the HGNC harmonisation "
                                 "of the cross-cohort analyses for a small number of genes.",
    "TCGA/TCGA.COAD.sampleMap_COAD_clinicalMatrix.txt": "UCSC Xena updates the TCGA clinical matrix occasionally. In "
                                 "addition, the copy used for the paper had been opened and re-saved in a spreadsheet "
                                 "program (dates and number formats rewritten, e.g. '26/02/2015' for '2015-2-26'). The "
                                 "fields the pipeline reads (sample_type, pathologic_stage, MSI status CDE_ID_3226963) "
                                 "were checked to be identical to the Xena file of 14 October 2026.",
}
KNOWN_DIFF = {  # recorded checksum not reproducible from the source; explained, but still an error unless allowed
    "TCGA/COAD_mc3_gene_level.txt": "The copy used for the paper had been opened and re-saved in a spreadsheet program "
                                    "(Windows line endings; gene symbols such as MARC1, SEP15 and FEB6 converted to dates). "
                                    "The pipeline uses only the number of mutated genes per sample, which was checked to "
                                    "be identical to the Xena file of 14 October 2026. Use --allow-mismatch to continue.",
}
GEO_ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc={acc}&form=text&view=full"


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def expected():
    out = {}
    for line in open(os.path.join(HERE, "PROVENANCE.md"), encoding="utf-8"):
        m = re.match(r"\| `data/raw/([^`]+)` \| [\d,]+ \| [^|]+ \| `([0-9a-f]{64})` \|", line)
        if m:
            out[m.group(1)] = m.group(2)
    return out


def fetch(url, dst, gunzip, tries=4):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    tmp = dst + ".part"
    for attempt in range(1, tries + 1):
        try:
            with requests.get(url, stream=True, timeout=(30, 120)) as r:   # connect / read timeouts (s)
                r.raise_for_status()
                with open(tmp, "wb") as fh:
                    for b in r.iter_content(1 << 20):
                        fh.write(b)
            break
        except requests.RequestException as e:
            print(f"          attempt {attempt} failed: {type(e).__name__}", flush=True)
            if attempt == tries:
                raise
    if gunzip:
        with gzip.open(tmp, "rb") as src, open(dst, "wb") as out:
            shutil.copyfileobj(src, out)
        os.remove(tmp)
    else:
        os.replace(tmp, dst)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--allow-mismatch", action="store_true", help="continue even if a non-volatile file differs")
    ap.add_argument("--check-only", action="store_true", help="do not download, only compare checksums")
    a = ap.parse_args()
    os.makedirs(RAW, exist_ok=True)
    if not a.check_only:
        for name, (url, gz) in FILES.items():
            dst = os.path.join(RAW, name)
            if os.path.exists(dst):
                print(f"present   {name}")
                continue
            print(f"download  {name}", flush=True)
            fetch(url, dst, gz)
        for gpl in PLATFORMS:   # the same GEO endpoint GEOparse uses (view=full), with timeouts and retries
            dst = os.path.join(RAW, f"{gpl}.txt")
            if os.path.exists(dst):
                print(f"present   {gpl}.txt")
                continue
            print(f"download  {gpl}.txt", flush=True)
            fetch(GEO_ACC.format(acc=gpl), dst, False)
    exp = expected()
    errors, warnings = [], []
    print("\nChecksums against data/PROVENANCE.md (files used for the paper):")
    for name in sorted(set(FILES) | {f"{g}.txt" for g in PLATFORMS} | {"GSE117606_cel_headers.csv"}):
        p = os.path.join(RAW, name)
        if not os.path.exists(p):
            print(f"  MISSING   {name}")
            errors.append(name)
            continue
        got = sha256(p)
        if exp.get(name) == got:
            print(f"  MATCH     {name}")
            continue
        detail = f"(expected {exp.get(name, 'not listed')[:12]}..., got {got[:12]}...)"
        if name in VOLATILE:
            print(f"  WARNING   {name} differs {detail}\n            {VOLATILE[name]}")
            warnings.append(name)
        else:
            print(f"  ERROR     {name} differs {detail}" + (f"\n            {KNOWN_DIFF[name]}" if name in KNOWN_DIFF else ""))
            errors.append(name)
    if warnings:
        print(f"\n{len(warnings)} file(s) from sources that are updated over time differ from the versions used for "
              "the paper (warning only). Results can differ slightly; the archived files give an exact reproduction "
              "(README, 'Data').")
    if errors:
        print(f"\n{len(errors)} other file(s) are missing or differ: {', '.join(errors)}.")
        if not a.allow_mismatch:
            sys.exit(1)
    if not warnings and not errors:
        print("\nAll files identical to the versions used for the paper.")


if __name__ == "__main__":
    main()
