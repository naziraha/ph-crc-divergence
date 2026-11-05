"""
Put every input file into data/raw/ and compare it with the SHA-256 recorded in data/PROVENANCE.md (the files used for
the paper).

Sources
  Archived copies (default for four files whose public versions change or differ; see data/PROVENANCE.md):
    Homo_sapiens.gene_info.gz, GPL570.txt, TCGA clinical matrix, TCGA MC3 gene-level mutations
    -> from the Zenodo deposit ([Zenodo DOI]); before it is published, pass --archive-dir <unpacked deposit folder>.
  Public sources (all other files):
    GEO series matrices GSE41258, GSE20916, GSE37364, GSE117606; GEO platform tables GPL96, GPL25373;
    UCSC Xena TCGA-COAD HiSeqV2 and GISTIC2 gene-level copy number.
  data/raw/GSE117606_cel_headers.csv is shipped with the repository and only checked here.

Checks
  - Every file is compared with data/PROVENANCE.md. With the archived copies, all files should MATCH, and any
    difference is an ERROR (exit code 1).
  - --from-source downloads the current public versions of the four archived files instead. They are known to differ
    (NCBI gene_info is regenerated daily; the GEO platform header lists newly submitted samples; the two TCGA files used
    for the paper had been re-saved in a spreadsheet program, without changing the fields used), so these differences
    are reported as WARNINGS with an explanation. Results can then differ slightly.
  - --allow-mismatch continues after an ERROR.

Usage: python data/download_data.py [--archive-dir DIR | --from-source] [--allow-mismatch] [--check-only]
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
ZENODO_FILES = ""   # base URL of the deposited files, e.g. "https://zenodo.org/records/<id>/files/"; set on publication
GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}nnn/{gse}/matrix/{gse}_series_matrix.txt.gz"
GEO_ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc={acc}&form=text&view=full"
XENA = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
SOURCES = {  # local name -> (public source URL, decompress gzip after download)
    **{f"{g}_series_matrix.txt.gz": (GEO_FTP.format(stub=g[:-3], gse=g), False)
       for g in ("GSE41258", "GSE20916", "GSE37364", "GSE117606")},
    **{f"{p}.txt": (GEO_ACC.format(acc=p), False) for p in ("GPL96", "GPL570", "GPL25373")},
    "Homo_sapiens.gene_info.gz": ("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz", False),
    "TCGA/TCGA.COAD.sampleMap_HiSeqV2.gz": (XENA + "TCGA.COAD.sampleMap%2FHiSeqV2.gz", False),
    "TCGA/TCGA.COAD.sampleMap_COAD_clinicalMatrix.txt": (XENA + "TCGA.COAD.sampleMap%2FCOAD_clinicalMatrix", False),
    "TCGA/COAD_mc3_gene_level.txt": (XENA + "mc3_gene_level%2FCOAD_mc3_gene_level.txt.gz", True),
    "TCGA/TCGA.COAD.sampleMap_Gistic2_all_data_by_genes.gz": (XENA + "TCGA.COAD.sampleMap%2FGistic2_CopyNumber_Gistic2_all_data_by_genes.gz", False),
}
ARCHIVED = ("Homo_sapiens.gene_info.gz", "GPL570.txt", "TCGA/TCGA.COAD.sampleMap_COAD_clinicalMatrix.txt",
            "TCGA/COAD_mc3_gene_level.txt")
_GPL = ("GEO appends every newly submitted sample and series that uses the platform to the header of the platform "
        "record, so the file changes over time. The pipeline reads only the probe annotation table (ID, Gene Symbol); "
        "for GPL570 it was checked to be identical to the version of 14 October 2026.")
CHANGES = {  # why the current public version of a file can differ from the one used for the paper
    "GPL96.txt": _GPL, "GPL570.txt": _GPL, "GPL25373.txt": _GPL,
    "Homo_sapiens.gene_info.gz": "NCBI regenerates gene_info daily (new genes, renamed symbols, retired IDs). A newer "
        "file can change the Entrez-to-symbol mapping of GSE117606 and the HGNC harmonisation of the cross-cohort "
        "analyses for a small number of genes.",
    "TCGA/TCGA.COAD.sampleMap_COAD_clinicalMatrix.txt": "The copy used for the paper had been re-saved in a spreadsheet "
        "program (dates and number formats rewritten); Xena may also update the file. The fields used (sample_type, "
        "pathologic_stage, MSI status CDE_ID_3226963) were checked to be identical to the Xena file of 14 October 2026.",
    "TCGA/COAD_mc3_gene_level.txt": "The copy used for the paper had been re-saved in a spreadsheet program (Windows line "
        "endings; gene symbols such as MARC1, SEP15 and FEB6 converted to dates). The only quantity used, the number of "
        "mutated genes per sample, was checked to be identical to the Xena file of 14 October 2026.",
}


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
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--archive-dir", help="folder with the unpacked Zenodo deposit (used instead of the Zenodo URL)")
    g.add_argument("--from-source", action="store_true", help="download current public versions of the archived files")
    ap.add_argument("--allow-mismatch", action="store_true", help="continue after a checksum ERROR")
    ap.add_argument("--check-only", action="store_true", help="do not download, only compare checksums")
    a = ap.parse_args()
    os.makedirs(RAW, exist_ok=True)
    if not a.check_only:
        for name, (url, gz) in SOURCES.items():
            dst = os.path.join(RAW, name)
            if os.path.exists(dst):
                print(f"present   {name}")
                continue
            if name in ARCHIVED and not a.from_source:
                if a.archive_dir:
                    print(f"archive   {name}  (from {a.archive_dir})", flush=True)
                    os.makedirs(os.path.dirname(dst), exist_ok=True)
                    src = os.path.join(a.archive_dir, name)
                    if not os.path.exists(src):          # the Zenodo deposit holds the files without sub-folders
                        src = os.path.join(a.archive_dir, os.path.basename(name))
                    shutil.copy2(src, dst)
                    continue
                if not ZENODO_FILES:
                    sys.exit(f"{name}: the archived copy is not yet published on Zenodo. Pass --archive-dir <folder with "
                             "the unpacked deposit>, or --from-source to use the current public version (results can "
                             "then differ slightly).")
                print(f"archive   {name}  (Zenodo)", flush=True)
                fetch(ZENODO_FILES + name.split("/")[-1] + "?download=1", dst, False)
                continue
            print(f"download  {name}", flush=True)
            fetch(url, dst, gz)
    exp = expected()
    errors, warnings = [], []
    print("\nChecksums against data/PROVENANCE.md (files used for the paper):")
    for name in sorted(set(SOURCES) | {"GSE117606_cel_headers.csv"}):
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
        known = (name not in ARCHIVED or a.from_source) and name in CHANGES
        if known:
            print(f"  WARNING   {name} differs {detail}\n            {CHANGES[name]}")
            warnings.append(name)
        else:
            print(f"  ERROR     {name} differs {detail}" + (f"\n            {CHANGES[name]}" if name in CHANGES else ""))
            errors.append(name)
    if warnings:
        print(f"\n{len(warnings)} file(s) differ from the versions used for the paper because their public sources "
              "change over time (warning only). Results can differ slightly; the archived copies give an exact "
              "reproduction.")
    if errors:
        print(f"\n{len(errors)} file(s) are missing or differ: {', '.join(errors)}.")
        if not a.allow_mismatch:
            sys.exit(1)
    if not warnings and not errors:
        print("\nAll files identical to the versions used for the paper.")


if __name__ == "__main__":
    main()
