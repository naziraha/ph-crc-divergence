"""
Addendum 6 section 7.1 - read CEL HEADER fields of GSE117606 (no intensities).

The GEO archive GSE117606_RAW.tar (path in the environment variable GSE117606_RAW_TAR; default
data/raw/GSE117606_RAW.tar) is read member by member;
for each Command Console ("Calvin") CEL.gz only the first 256 KB of decompressed data are read and
the generic data header (with its parent headers, e.g. the scan acquisition) is parsed. The data
groups / intensity block are never read. Output: data/raw/GSE117606_cel_headers.csv.
"""
import gzip
import io
import os
import struct
import tarfile

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAR = os.environ.get("GSE117606_RAW_TAR", os.path.join(ROOT, "data", "raw", "GSE117606_RAW.tar"))
OUT = os.path.join(ROOT, "data", "raw", "GSE117606_cel_headers.csv")
KEEP = ("date", "barcode", "plate", "scanner", "well", "position", "array-type", "time")


class Reader:
    def __init__(self, b):
        self.b, self.p = b, 0

    def i32(self):
        v = struct.unpack(">i", self.b[self.p:self.p + 4])[0]; self.p += 4; return v

    def u8(self):
        v = self.b[self.p]; self.p += 1; return v

    def astr(self):
        n = self.i32(); s = self.b[self.p:self.p + n]; self.p += n; return s.decode("latin-1")

    def wstr(self):
        n = self.i32(); s = self.b[self.p:self.p + 2 * n]; self.p += 2 * n
        return s.decode("utf-16-be", "ignore").rstrip("\x00")

    def raw(self):
        n = self.i32(); s = self.b[self.p:self.p + n]; self.p += n; return s


def value(raw, typ):
    if typ.startswith("text/plain"):
        return raw.decode("utf-16-be", "ignore").rstrip("\x00").strip()
    if typ.startswith("text/ascii"):
        return raw.decode("latin-1", "ignore").rstrip("\x00").strip()
    if "integer-32" in typ and len(raw) >= 4:
        return struct.unpack(">i", raw[:4])[0]
    if "float" in typ and len(raw) >= 4:
        return struct.unpack(">f", raw[:4])[0]
    return raw.hex()[:32]


def header(r, depth, out):
    dtype = r.astr(); r.astr()                 # data type id, unique file id
    out.setdefault(f"h{depth}_type", dtype)
    out.setdefault(f"h{depth}_datetime", r.wstr())
    r.wstr()                                    # locale
    for _ in range(r.i32()):
        name, raw, typ = r.wstr(), r.raw(), r.wstr()
        if any(k in name.lower() for k in KEEP):
            out.setdefault(name, value(raw, typ))
    for _ in range(r.i32()):                    # parent headers (e.g. scan acquisition)
        header(r, depth + 1, out)


rows = []
with tarfile.open(TAR) as t:
    for m in t.getmembers():
        if not m.name.endswith(".CEL.gz"):
            continue
        comp = t.extractfile(m).read(300_000)            # compressed bytes of the header region only
        data = gzip.GzipFile(fileobj=io.BytesIO(comp)).read(262_144)
        r = Reader(data)
        assert r.u8() == 59, m.name                     # Calvin magic
        r.u8(); r.i32(); r.i32()                        # version, n data groups, first group position
        rec = {"file": m.name, "gsm": m.name.split("_")[0]}
        header(r, 0, rec)
        rows.append(rec)
df = pd.DataFrame(rows)
cols = ["gsm", "file", "affymetrix-array-type", "affymetrix-scan-date", "affymetrix-plate-barcode",
        "affymetrix-ArrayPlateGUID", "affymetrix-PlateScanGUID", "affymetrix-plate-peg-wellposition",
        "affymetrix-array-barcode", "affymetrix-Hyb-Start-Time", "affymetrix-Hyb-Stop-Time",
        "affymetrix-Fluidics-Date", "affymetrix-scanner-id", "affymetrix-scanner-serialnumber"]
df = df[[c for c in cols if c in df.columns]]
df.to_csv(OUT, index=False)
print(f"{len(df)} CEL headers -> {OUT}")
print("columns:", list(df.columns))
print(df.head(3).T.to_string())
