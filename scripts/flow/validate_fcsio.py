# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: One-off check of fcsio.py against the fcsparser package.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""One-off validation of fcsio.py against fcsparser (decision 2.9): events and TEXT identical on 3 tubes.

Run with the kmer_analysis python (has fcsparser), isolated and without bytecode writes:
  ~/miniconda3/envs/kmer_analysis/bin/python -I -B scripts/flow/validate_fcsio.py <fcs> <fcs> <fcs>
"""
import os
import sys

import fcsparser
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fcsio  # noqa: E402

ok_all = True
for f in sys.argv[1:]:
    meta, df = fcsparser.parse(f, reformat_meta=False)
    names, X, kv = fcsio.read(f)
    same_shape = X.shape == df.values.shape
    same_vals = same_shape and np.array_equal(X, df.values.astype(float))
    same_names = list(df.columns) == names
    mk = {k.upper(): str(v) for k, v in meta.items() if not k.startswith("__")}
    keys = [k for k in mk if k in kv]
    diff = [k for k in keys if mk[k].strip() != kv[k].strip()]
    missing = [k for k in mk if k not in kv]
    ok = same_vals and same_names and not diff and not missing
    ok_all &= ok
    print(f"{os.path.basename(f)}: events {X.shape} values identical {same_vals}; names identical {same_names}; "
          f"TEXT keys compared {len(keys)}, differing {diff}, missing in fcsio {missing} -> {'OK' if ok else 'MISMATCH'}")
print("ALL IDENTICAL" if ok_all else "MISMATCH FOUND")
