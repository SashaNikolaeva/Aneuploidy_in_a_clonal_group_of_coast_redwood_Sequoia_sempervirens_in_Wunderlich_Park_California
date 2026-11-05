#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Freezes the written plan of the analysis run (prereg/finalize_PREREG.tsv) before any data
# are opened (Methods: the tests and thresholds were written down before the clone statistics were computed).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""Step 0: freeze the PREREG table of the 2026-10-26 decision as PREREG.tsv (read-only).

Copies the markdown table between <!-- PREREG:BEGIN --> / <!-- PREREG:END --> verbatim
(cell text unchanged; only the markdown pipes become tabs), asserts the row ids,
writes PREREG.sha256 (and the decision's sha256), chmod 0444, refuses to overwrite.
Run BEFORE any input data file is opened.
"""
import hashlib, os, sys, datetime

REPO = "/space/s1/sashanikolaeva/Scripts/GitHub/redwood-clone-paper"
DEC = f"{REPO}/analysis/2026-10-26_finalize_analyses_decision.md"
OUT = f"{REPO}/analysis/finalize_run"
EXPECTED_DEC_SHA = "f24c1d5b26fcc81ea3b291b15f98b979fe9fe691426f16021519a2a29a3a2b34"
IDS = ["R0g", "R1a", "R1b", "R1c-g", "R1c", "R2", "R2-WPR", "R2-364", "R3", "R4", "R5", "R6-oct", "R6-run"]

def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()

def main():
    tsv, shaf = f"{OUT}/PREREG.tsv", f"{OUT}/PREREG.sha256"
    if os.path.exists(tsv) or os.path.exists(shaf):
        sys.exit("REFUSE: PREREG.tsv / PREREG.sha256 already exist; changes go to DEVIATIONS.tsv")
    d = sha(DEC)
    if d != EXPECTED_DEC_SHA:
        sys.exit(f"STOP: decision sha256 {d} != recorded {EXPECTED_DEC_SHA}")
    lines = open(DEC).read().splitlines()
    b = lines.index("<!-- PREREG:BEGIN -->"); e = lines.index("<!-- PREREG:END -->")
    rows = [l for l in lines[b + 1:e] if l.strip()]
    out = []
    for i, l in enumerate(rows):
        cells = [c.strip() for c in l.strip().strip("|").split("|")]
        if i == 1:  # markdown separator row
            assert all(set(c) <= set("-:") for c in cells); continue
        out.append(cells)
    hdr = out[0]
    assert hdr[0] == "id" and hdr[-1] == "already seen?" and len(hdr) == 12, hdr
    for r in out[1:]:
        assert len(r) == 12, r
        assert "\t" not in "".join(r)
    ids = [r[0] for r in out[1:]]
    assert ids == IDS, ids
    with open(tsv, "w") as f:
        for r in out:
            f.write("\t".join(r) + "\n")
    h = sha(tsv)
    with open(shaf, "w") as f:
        f.write(f"{h}  PREREG.tsv\n{d}  ../2026-10-26_finalize_analyses_decision.md\n")
        f.write(f"# frozen {datetime.datetime.now().isoformat(timespec='seconds')} before any input data file was opened\n")
    os.chmod(tsv, 0o444); os.chmod(shaf, 0o444)
    print(f"PREREG.tsv frozen: {len(ids)} rows, sha256 {h}")

if __name__ == "__main__":
    main()
