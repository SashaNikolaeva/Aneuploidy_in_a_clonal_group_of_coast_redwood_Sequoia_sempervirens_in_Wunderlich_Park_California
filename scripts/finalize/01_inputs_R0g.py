# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Records input checksums and the fixed sample sets, and checks that the inputs reproduce the
# earlier chromosome-4 value before the run continues.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""Step 1: input sha256, fixed sample-set asserts (decision 2.0), and the R0g equivalence gate.

R0g (PREREG): clone chr4 genet c_hat = 5.024 +- 0.005, z vs PANEL_ALL = -14.4 +- 0.1, n(PANEL_ALL) = 223.
  c_hat = 6 * R / median(PANEL_ALL R); genet value = mean of the 9 ramets' R; z = (genet - mean_P) / sd_P.
On failure: write STOP.md and exit 2 (stop rule S0).
Outputs: input_hashes.tsv, sample_sets.tsv, R0g_gate.tsv, 01_inputs_R0g.log (stdout).
"""
import os
import statistics as st
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *  # noqa

ensure_dirs()
hf = f"{OUT}/input_hashes.tsv"
rows = []
for k, p in INPUTS.items():
    if not os.path.exists(p):
        open(f"{OUT}/STOP.md", "w").write(f"# STOP (S1)\n\nInput missing: {k} = {p}\n")
        sys.exit(f"STOP S1: missing input {p}")
    rows.append((k, p, sha256(p), os.path.getsize(p)))
if os.path.exists(hf):   # re-run: hashes must not have changed (S1)
    old = {l.split("\t")[0]: l.split("\t")[2] for l in open(hf).read().splitlines()[1:]}
    for k, p, h, _ in rows:
        if k in old and old[k] != h:
            open(f"{OUT}/STOP.md", "w").write(f"# STOP (S1)\n\nInput hash changed: {k}\n")
            sys.exit(f"STOP S1: hash changed for {k}")
with open(hf, "w") as f:
    f.write("key\tpath\tsha256\tbytes\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")
print(f"input hashes written: {len(rows)} files")

S = sample_sets()
rc = S["rc"]
print(f"SEGI counts rows: {len(rc)}; > 1 M reads: {len(S['big'])}")
missing = [s for s in CLONE if s not in rc]
n = {k: len(S[k]) for k in ("ref229", "PANEL_ALL", "PANEL_LOCAL", "PANEL_SOUTH")}
print(f"CLONE 9 present: {not missing}; ref229 {n['ref229']}; PANEL_ALL {n['PANEL_ALL']}; LOCAL {n['PANEL_LOCAL']}; "
      f"SOUTH {n['PANEL_SOUTH']} (reported, expected ~43); centroid {S['centroid'][0]:.5f}, {S['centroid'][1]:.5f}")
inv = S["inv"]
gids = [inv[s]["genet_id"] for s in S["PANEL_ALL"]]
dups = sorted({g for g in gids if gids.count(g) > 1})
print(f"PANEL_ALL genet ids unique: {not dups} {dups if dups else ''}")
assert all(inv[s]["genet_id"] == "G000" for s in CLONE)
tech_ok = all(a in rc and b in rc for a, b in TECH)
gain_ok = all(t in S["PANEL_ALL"] for t in GAINCTRL)
print(f"TECH pairs present: {tech_ok}; GAINCTRL in PANEL_ALL: {gain_ok}; WPR in panel {WPR in S['PANEL_ALL']}; "
      f"C1 in panel {C1 in S['PANEL_ALL']}; 14_364 in panel {S364 in S['PANEL_ALL']}")
fail = []
if missing: fail.append(f"clone ramets missing {missing}")
N_EXP = 223 if SET == "prereg" else 224   # popgen224 (DV14): exactly the popgen 224
if n["ref229"] != 229: fail.append(f"ref229 = {n['ref229']}")
if n["PANEL_ALL"] != N_EXP: fail.append(f"PANEL_ALL = {n['PANEL_ALL']}")
if SET == "prereg" and n["PANEL_LOCAL"] != 35: fail.append(f"PANEL_LOCAL = {n['PANEL_LOCAL']}")
print(f"set {SET}: PANEL_ALL {n['PANEL_ALL']} (expected {N_EXP}); PANEL_LOCAL {n['PANEL_LOCAL']}; PANEL_SOUTH {n['PANEL_SOUTH']}")
if dups: fail.append(f"duplicate genets in PANEL_ALL {dups}")
if not tech_ok or not gain_ok: fail.append("TECH/GAINCTRL")
with open(f"{OUT}/sample_sets.tsv", "w") as f:
    f.write("set\tsample\tgenet_id\tlat\tlon\tregion_lat\tpcr_pool\tpcr_group\n")
    for k in ("PANEL_ALL", "PANEL_LOCAL", "PANEL_SOUTH"):
        for s in S[k]:
            r = inv.get(s, {})
            f.write(f"{k}\t{s}\t{r.get('genet_id','')}\t{r.get('lat','')}\t{r.get('lon','')}\t{r.get('region_lat','')}"
                    f"\t{r.get('pcr_pool','')}\t{r.get('pcr_group','')}\n")
    for s in CLONE:
        r = inv[s]
        f.write(f"CLONE\t{s}\t{r['genet_id']}\t{r['lat']}\t{r['lon']}\t{r['region_lat']}\t{r['pcr_pool']}\t{r['pcr_group']}\n")
    for a, b in TECH:
        f.write(f"TECH\t{a}|{b}\t\t\t\t\t\t\n")
    for t, c in GAINCTRL.items():
        f.write(f"GAINCTRL\t{t}:{c}\t\t\t\t\t\t\n")
if fail:
    open(f"{OUT}/STOP.md", "w").write("# STOP (S1): sample-set mismatch\n\n" + "\n".join(f"- {x}" for x in fail) + "\n")
    sys.exit("STOP S1: " + "; ".join(fail))

# ---- R0g
P = S["PANEL_ALL"]
R = {s: segi_ratio(rc, s, "chr4") for s in S["big"]}
pv = [R[s] for s in P]
med, mu, sd = st.median(pv), st.mean(pv), st.stdev(pv)
g = st.mean(R[s] for s in CLONE)
chat = 6 * g / med
z = (g - mu) / sd
t, df, pl, _ = pred_t(g, pv)
ok = abs(chat - 5.024) <= 0.005 and abs(z - (-14.4)) <= 0.1 and len(P) == N_EXP
if SET == "popgen224":   # the R0g targets belong to the 223-tree prereg panel: logged, not fatal (DV14)
    print("popgen224: R0g targets are for the prereg panel; value logged only")
print(f"\nR0g: genet c_hat {chat:.4f} (target 5.024 +- 0.005); z {z:+.3f} (target -14.4 +- 0.1); n {len(P)} (223); "
      f"[prediction t {t:+.2f}, df {df}] -> {'PASS' if ok else 'FAIL'}")
with open(f"{OUT}/R0g_gate.tsv", "w") as f:
    f.write("statistic\tvalue\ttarget\ttolerance\tpass\n")
    f.write(f"c_hat_genet_chr4\t{chat:.4f}\t5.024\t0.005\t{abs(chat-5.024)<=0.005}\n")
    f.write(f"z_vs_PANEL_ALL\t{z:.3f}\t-14.4\t0.1\t{abs(z+14.4)<=0.1}\n")
    f.write(f"n_PANEL_ALL\t{len(P)}\t223\t0\t{len(P)==223}\n")
    f.write(f"panel_median_R\t{med:.5f}\t\t\t\n")
    f.write(f"panel_mean_R\t{mu:.5f}\t\t\t\n")
    f.write(f"panel_sd_R\t{sd:.5f}\t\t\t\n")
    f.write(f"GATE\t{'PASS' if ok else 'FAIL'}\t\t\t{ok}\n")
if not ok and SET == "prereg":
    open(f"{OUT}/STOP.md", "w").write(
        "# STOP (S0): R0g failed - input drift: SEGI counts do not reproduce rev01\n\n"
        f"c_hat {chat:.4f} (5.024 +- 0.005), z {z:+.3f} (-14.4 +- 0.1), n {len(P)} (223)\n")
    sys.exit(2)
