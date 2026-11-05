# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Freezes the parameters of a flow-cytometry run before the pipeline is run on it.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""Freeze PREREG_flow_<run>.tsv (step 8) BEFORE the pipeline is run on that run, and in a designed run BEFORE the
blinded tube key is unblinded. Writes parameters (defaults of flowpipe.DEFAULTS, overridable with key=value), the
pre-registered outcomes (proposal 4c) and the genet-level power table (review flow_04), then the sha256 and chmod 0444.
Refuses to overwrite: any later change goes into the run's DEVIATIONS.

  python3 -I freeze_prereg.py <run> <out_dir> [--purpose "..."] [key=value ...]
"""
import argparse
import datetime
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flowpipe  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run"); ap.add_argument("out_dir")
    ap.add_argument("--purpose", default="designed run: clone vs euploid delta within stratum")
    ap.add_argument("set", nargs="*")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    p = os.path.join(a.out_dir, f"PREREG_flow_{a.run}.tsv")
    if os.path.exists(p):
        sys.exit(f"REFUSE: {p} exists and is frozen")
    P = dict(flowpipe.DEFAULTS)
    for kv in a.set:
        k, v = kv.split("=", 1)
        if k not in P:
            sys.exit(f"unknown parameter {k}")
        P[k] = v
    with open(p, "w") as f:
        f.write("section\tkey\tvalue\n")
        f.write(f"run\trun\t{a.run}\nrun\tpurpose\t{a.purpose}\n")
        for k, v in P.items():
            f.write(f"param\t{k}\t{v}\n")
        f.write("outcome\tunit\tindependent chop (tubes of one chop averaged); tree = cluster; contrast within run x voltage stratum\n")
        f.write("outcome\tstatistic\tdelta % = (clone chop mean - mean of euploid tree means) / mean of euploid tree means; "
                "95% CI by two-level bootstrap (trees, then chops); PASS-B and above only\n")
        f.write(f"outcome\tconsistent\tCI excludes 0 and point in [{P['out_consistent_lo']}, {P['out_consistent_hi']}]% "
                "-> 'consistent with one-copy loss' (not chr4-specific)\n")
        f.write(f"outcome\tdoes_not_support\tCI includes 0 with half-width <= {P['out_null_halfwidth']}%\n")
        f.write(f"outcome\tinconsistent\tdelta < {P['out_inconsistent_below']}% with CI excluding {P['expected_delta_pct']}% "
                "-> artefact or larger loss; re-examine\n")
        f.write("outcome\tinconclusive\tanything else\n")
        f.write("outcome\tdoes_not_license\twhich chromosome; a copy number\n")
        for k, need, smax in flowpipe.power_table():
            f.write(f"power\tk_euploid_trees={k}\trequired delta/sigma_T {need:.2f}; max sigma_T for delta 1.54% = {smax:.2f}%\n")
        f.write(f"frozen\tat\t{datetime.datetime.now().isoformat(timespec='seconds')}\n")
    h = hashlib.sha256(open(p, "rb").read()).hexdigest()
    with open(p + ".sha256", "w") as f:
        f.write(f"{h}  {os.path.basename(p)}\n")
    os.chmod(p, 0o444); os.chmod(p + ".sha256", 0o444)
    print(f"frozen {p} sha256 {h}")


if __name__ == "__main__":
    main()
