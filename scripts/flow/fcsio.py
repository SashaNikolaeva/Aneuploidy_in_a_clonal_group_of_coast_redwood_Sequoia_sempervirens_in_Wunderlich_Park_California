# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Minimal FCS file reader used by the flow-cytometry pipeline.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""Minimal FCS 2.0 / 3.0 / 3.1 reader (numpy only; no fcsparser). Promoted from analysis/finalize_review/fcsio.py.

text(path)  -> dict of TEXT keywords (keys upper-cased; '$P1V' etc.), plus '_HDR_DATA' = (start, end) from the header
read(path)  -> (names, events x parameters float array, kv)
Supports $DATATYPE F / D / I (16- or 32-bit, uniform width), either byte order, and a doubled delimiter inside values.
Validated against fcsparser 0.2.8 (kmer_analysis env) by validate_fcsio.py.
"""
import numpy as np


def text(path):
    with open(path, "rb") as fh:
        h = fh.read(58)
        if not h[:3] == b"FCS":
            raise ValueError(f"{path}: not an FCS file")
        ts, te = int(h[10:18]), int(h[18:26])
        ds = int(h[26:34].strip() or 0)
        de = int(h[34:42].strip() or 0)
        fh.seek(ts)
        raw = fh.read(te - ts + 1).decode("latin-1")
    d = raw[0]
    body = raw[1:]
    if body.endswith(d):
        body = body[:-1]
    # a doubled delimiter is an escaped delimiter inside a value
    tok = body.replace(d + d, "\x00").split(d)
    tok = [t.replace("\x00", d) for t in tok]
    kv = {}
    for i in range(0, len(tok) - 1, 2):
        kv[tok[i].strip().upper()] = tok[i + 1]
    kv["_HDR_DATA"] = (ds, de)
    return kv


def read(path):
    kv = text(path)
    n = int(kv["$PAR"])
    tot = int(kv["$TOT"].strip())
    ds, de = kv["_HDR_DATA"]
    if ds == 0:
        ds, de = int(kv["$BEGINDATA"]), int(kv["$ENDDATA"])
    dt = kv["$DATATYPE"].strip().upper()
    bo = kv["$BYTEORD"].strip()
    bits = {int(kv[f"$P{i}B"]) for i in range(1, n + 1)}
    if len(bits) != 1:
        raise ValueError(f"{path}: mixed parameter bit widths {bits} not supported")
    b = bits.pop()
    endian = "<" if bo.startswith("1,2") else ">"
    code = {"F": "f4", "D": "f8", "I": {16: "u2", 32: "u4"}.get(b)}[dt]
    with open(path, "rb") as fh:
        fh.seek(ds)
        buf = fh.read(de - ds + 1)[: tot * n * b // 8]
    arr = np.frombuffer(buf, dtype=endian + code).reshape(tot, n).astype(float)
    names = [kv[f"$P{i}N"].strip() for i in range(1, n + 1)]
    return names, arr, kv


def param_meta(kv):
    """Per-parameter name -> dict(V=voltage str, R=range float, G=gain str)."""
    n = int(kv["$PAR"])
    out = {}
    for i in range(1, n + 1):
        out[kv[f"$P{i}N"].strip()] = dict(V=kv.get(f"$P{i}V", "").strip(), R=float(kv.get(f"$P{i}R", "nan")),
                                         G=kv.get(f"$P{i}G", "").strip())
    return out
