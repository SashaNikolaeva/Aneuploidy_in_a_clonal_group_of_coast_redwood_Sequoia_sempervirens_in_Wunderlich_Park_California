# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Builds the sample sheet of a flow-cytometry run from its FCS files.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""Build a skeleton sample sheet from a folder of FCS files, for the lab / Sasha to complete (step 0).

  python3 -I make_sample_sheet.py <data_dir> <run_name> <out.tsv> [--standard SEGI] [--scan '<glob>' ...]

One row per unique FCS file (byte-identical copies, e.g. a re-upload in a sub-folder, are listed once; the copy is
named in notes). Filled automatically: file (relative to data_dir), tube, run_date ($DATE), voltage_set, role
(standard_only = the tube name is the standard code alone; sample = contains '+' or no standard at all).
Left for the lab: tree (WP ramet ID / WPR / C1 / other tree ID), tree_class (clone / euploid / unknown), tissue,
buffer, chop_id (one ID per independent chop; re-runs of one chop share it), rerun_of, blinded_code.
Unfilled values are UNKNOWN; run_flow.py never computes a clone-vs-euploid contrast while any tree is UNKNOWN.
"""
import argparse
import glob
import hashlib
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fcsio  # noqa: E402

COLS = ["file", "tube", "run", "run_date", "role", "tree", "tree_class", "tissue", "buffer", "standard", "chop_id",
        "rerun_of", "voltage_set", "blinded_code", "dna_channel", "include", "notes"]


def tube_name(path, kv):
    t = kv.get("TUBE NAME", "").strip()
    if not t:
        t = re.sub(r"^Specimen_\d+_", "", os.path.basename(path)[:-4])
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data_dir"); ap.add_argument("run"); ap.add_argument("out")
    ap.add_argument("--standard", default="SEGI", help="code of the in-tube standard in tube names ('' = none)")
    ap.add_argument("--scan", action="append", default=None, help="glob(s) relative to data_dir (default '**/*.fcs')")
    a = ap.parse_args()
    if os.path.exists(a.out):
        sys.exit(f"REFUSE: {a.out} exists (edit it instead)")
    files = sorted({f for g in (a.scan or ["**/*.fcs"]) for f in glob.glob(os.path.join(a.data_dir, g), recursive=True)})
    seen, rows = {}, []
    for f in sorted(files, key=lambda p: (p.count(os.sep), p)):   # top-level copies first
        md5 = hashlib.md5(open(f, "rb").read()).hexdigest()
        rel = os.path.relpath(f, a.data_dir)
        if md5 in seen:
            seen[md5]["notes"] += f" identical copy: {rel};"
            continue
        kv = fcsio.text(f)
        names = [kv[f"$P{i}N"].strip() for i in range(1, int(kv["$PAR"]) + 1)]
        pm = fcsio.param_meta(kv)
        dna = next((c for c in ["BUV 496-A", "Pacific Blue-A", "BV421-A"] if c in names), "")
        t = tube_name(f, kv)
        std = a.standard
        if std and re.fullmatch(rf"{re.escape(std)}(_\d+)?", re.sub(r"_\d{3}$", "", t)):
            role = "standard_only"
        elif std and std in t:
            role = "sample"
        else:
            role, std = "sample", "none"
        vs = f"FSC{pm.get('FSC-A', {}).get('V', '')}/SSC{pm.get('SSC-A', {}).get('V', '')}/DNA{pm.get(dna, {}).get('V', '')}"
        r = dict(file=rel, tube=t, run=a.run, run_date=kv.get("$DATE", "").strip(), role=role,
                 tree="standard" if role == "standard_only" else "UNKNOWN",
                 tree_class="standard" if role == "standard_only" else "unknown", tissue="UNKNOWN", buffer="UNKNOWN",
                 standard=std, chop_id="UNKNOWN", rerun_of="", voltage_set=vs, blinded_code="", dna_channel="",
                 include="yes", notes="AUTO-GENERATED skeleton;")
        seen[md5] = r
        rows.append(r)
    with open(a.out, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in sorted(rows, key=lambda r: (r["run_date"], r["file"])):
            fh.write("\t".join(str(r[c]) for c in COLS) + "\n")
    print(f"{len(rows)} unique tubes ({len(files)} files) -> {a.out}")


if __name__ == "__main__":
    main()
