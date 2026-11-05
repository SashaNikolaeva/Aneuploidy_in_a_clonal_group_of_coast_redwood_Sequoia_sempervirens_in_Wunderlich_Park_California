# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Runs the flow-cytometry pipeline on one run folder.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""Run the flow genome-size pipeline on one run folder (see README.md).

  python3 -I run_flow.py --data <data_dir> --sheet <sample_sheet.tsv> --run <run> --out <out_dir> --figdir <fig_dir>
                         [--scan '<glob relative to data_dir>' ...]

Requires <out_dir>/PREREG_flow_<run>.tsv frozen by freeze_prereg.py (parameters are read from it).
Refuses FCS files in <data_dir> (or in the --scan globs) that are not in the sheet (byte-identical copies of listed files are allowed).
Writes to <out_dir>/<run>/: ingest.tsv, standards.tsv, per_tube.tsv, per_stratum.tsv, per_tree.tsv, contrast.tsv,
run.log; figures to <fig_dir>/<run>/.
"""
import argparse
import csv
import glob
import hashlib
import math
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flowpipe as fp  # noqa: E402

INK, INK2, GRID = "#0b0b0b", "#52514e", "#d6d5d0"
TIERCOL = {"PASS-A": "#2a78d6", "PASS-B": "#1baf7a", "FAIL": "#9a9994"}
CLASSMK = {"clone": "o", "euploid": "s", "unknown": "D", "standard": "^"}


class Log:
    def __init__(self, p):
        self.f = open(p, "w")

    def __call__(self, *a):
        s = " ".join(str(x) for x in a)
        print(s)
        self.f.write(s + "\n")


def read_prereg(p):
    P = dict(fp.DEFAULTS)
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["section"] == "param":
                P[r["key"]] = fp.num(r["value"])
    return P


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def fmt(v):
    if isinstance(v, float):
        return "NA" if math.isnan(v) else f"{v:.4f}"
    return str(v)


def main():
    ap = argparse.ArgumentParser()
    for k in ("data", "sheet", "run", "out", "figdir"):
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--scan", action="append", default=None,
                    help="glob(s) relative to --data that belong to this run (default '**/*.fcs'); repeatable")
    a = ap.parse_args()
    pre = os.path.join(a.out, f"PREREG_flow_{a.run}.tsv")
    if not os.path.exists(pre) or os.access(pre, os.W_OK):
        sys.exit(f"REFUSE: {pre} missing or not frozen (run freeze_prereg.py first)")
    P = read_prereg(pre)
    od, fd = os.path.join(a.out, a.run), os.path.join(a.figdir, a.run)
    os.makedirs(od, exist_ok=True); os.makedirs(fd, exist_ok=True)
    L = Log(os.path.join(od, "run.log"))
    L(f"run {a.run}; data {a.data}; sheet {a.sheet} (sha256 {hashlib.sha256(open(a.sheet,'rb').read()).hexdigest()}); "
      f"PREREG {pre} (sha256 {hashlib.sha256(open(pre,'rb').read()).hexdigest()})")
    sheet = [r for r in csv.DictReader(open(a.sheet), delimiter="\t") if r["include"].strip().lower() == "yes"]
    listed = {os.path.normpath(r["file"]) for r in sheet}
    lmd5 = {md5(os.path.join(a.data, r["file"])) for r in sheet}
    allf = sorted({os.path.relpath(f, a.data) for g in (a.scan or ["**/*.fcs"])
                   for f in glob.glob(os.path.join(a.data, g), recursive=True)})
    unlisted = [f for f in allf if os.path.normpath(f) not in listed and md5(os.path.join(a.data, f)) not in lmd5]
    if unlisted:
        sys.exit(f"REFUSE: FCS files not in the sample sheet: {unlisted}")
    L(f"tubes in sheet (include=yes): {len(sheet)}; FCS files in folder: {len(allf)} (identical copies allowed)")

    # ---- 1 ingest
    tubes = []
    for r in sheet:
        meta, A, H, ceil = fp.ingest(os.path.join(a.data, r["file"]), r.get("dna_channel", "").strip())
        tubes.append(dict(sheet=r, meta=meta, A=A, H=H, ceil=ceil,
                          stratum=(meta["date"], meta["dna_channel"], meta["voltage_set"])))
    with open(os.path.join(od, "ingest.tsv"), "w") as f:
        keys = list(tubes[0]["meta"])
        f.write("file\ttube\trole\t" + "\t".join(keys) + "\tsheet_voltage_set_matches\n")
        for t in tubes:
            m = t["meta"]
            match = (t["sheet"]["voltage_set"] in ("", m["voltage_set"]))
            f.write(f"{t['sheet']['file']}\t{t['sheet']['tube']}\t{t['sheet']['role']}\t" + "\t".join(fmt(m[k]) for k in keys)
                    + f"\t{match}\n")
    strata = sorted({t["stratum"] for t in tubes})
    L(f"strata (date x DNA channel x voltage set): {len(strata)}")
    for s in strata:
        L(f"  {s}: {sum(t['stratum'] == s for t in tubes)} tubes")

    # ---- 2 standards / gates
    refs = {}
    with open(os.path.join(od, "standards.tsv"), "w") as f:
        f.write("stratum\ttube\tG1\tCV\tn\tAH_lo\tAH_hi\n")
        for s in strata:
            st_ = [t for t in tubes if t["stratum"] == s and t["sheet"]["role"] == "standard_only"]
            got = []
            for t in st_:
                r = fp.std_reference(t["A"], t["H"], t["ceil"], P)
                if r:
                    got.append(r)
                    f.write(f"{'|'.join(s)}\t{t['sheet']['tube']}\t{r['g1']:.1f}\t{r['cv']:.2f}\t{r['n']}\t"
                            + (f"{r['ah_window'][0]:.4f}\t{r['ah_window'][1]:.4f}" if r["ah_window"] else "NA\tNA") + "\n")
                    L(f"standard-only {t['sheet']['tube']} [{s[0]} {s[2]}]: G1 {r['g1']:.0f} CV {r['cv']:.1f}% n {r['n']}")
            if got:
                win = None
                if all(g["ah_window"] for g in got):
                    win = (min(g["ah_window"][0] for g in got), max(g["ah_window"][1] for g in got))
                refs[s] = dict(g1=float(np.mean([g["g1"] for g in got])), ah=win, fallback=False,
                               n_std=len(got), src="standard-only tube(s)")
    for s in strata:
        if s in refs:
            continue
        std_codes = {t["sheet"]["standard"] for t in tubes if t["stratum"] == s}
        if std_codes == {"none"}:
            refs[s] = None
            L(f"stratum {s}: no standard in any tube -> ingest and ceiling QC only")
            continue
        cand = [(k, v) for k, v in refs.items() if v and k[1] == s[1] and not v["fallback"]]
        if not cand:
            refs[s] = None
            L(f"stratum {s}: no standard-only tube in the run for channel {s[1]} -> cannot gate")
            continue
        cand.sort(key=lambda kv: (kv[0][0] != s[0], kv[0]))
        k, v = cand[0]
        V, Vr = float(s[2].split("DNA")[1]), float(k[2].split("DNA")[1])
        refs[s] = dict(g1=v["g1"] * (V / Vr) ** P["gain_exp"], ah=None, fallback=True, n_std=0,
                       src=f"scaled from {k[0]} {k[2]} by (V/{Vr:.0f})^{P['gain_exp']}")
        L(f"stratum {s}: no standard-only tube; gate_fallback reference G1 {refs[s]['g1']:.0f} ({refs[s]['src']})")

    # ---- 3-5 per tube
    res = []
    for t in tubes:
        if t["sheet"]["role"] != "sample":
            continue
        ref = refs[t["stratum"]]
        row = dict(file=t["sheet"]["file"], tube=t["sheet"]["tube"], date=t["meta"]["date"], voltage_set=t["meta"]["voltage_set"],
                   dna_channel=t["meta"]["dna_channel"], events=t["meta"]["events"], tree=t["sheet"]["tree"],
                   tree_class=t["sheet"]["tree_class"], chop_id=t["sheet"]["chop_id"], ceil_dna_pct=t["meta"]["ceil_dna_pct"])
        if ref is None:
            row.update(tier="NO-STANDARD", tier_reason="no in-tube standard / no reference", flags="ingest_only",
                       ref_g1=np.nan)
            res.append((row, None, t)); continue
        r = fp.analyse_tube(t["A"], t["H"], t["ceil"], ref["g1"], ref["ah"], P)
        if not ref["fallback"] and not math.isnan(r[f"std_g1_{r['ratio_window_used'] or 'w1'}"]):
            r["suppression"] = r[f"std_g1_{r['ratio_window_used'] or 'w1'}"] / ref["g1"]
        if ref["fallback"]:
            r["flags"].append("gate_fallback")
        tr, why = fp.tier(r, P, suppression_known=not ref["fallback"])
        row.update({k: v for k, v in r.items() if k not in ("x", "fits", "flags")})
        row.update(tier=tr, tier_reason=why, flags=";".join(r["flags"]), ref_g1=ref["g1"])
        row["twoC_pg_assembly_std"] = row["ratio"] * P["std_2C_pg_assembly"] if not math.isnan(row["ratio"]) else np.nan
        lit = P["std_2C_pg_literature"]
        row["twoC_pg_literature_std"] = row["ratio"] * lit if isinstance(lit, float) and not math.isnan(row["ratio"]) else "NA"
        res.append((row, r, t))
    cols = ["file", "tube", "date", "voltage_set", "dna_channel", "events", "tree", "tree_class", "chop_id", "tier",
            "tier_reason", "flags", "ref_g1", "gated_events", "singlets", "std_g1_w1", "std_cv_w1", "std_n_w1", "red_g1_w1",
            "red_cv_w1", "red_n_w1", "std_g1_w2", "std_cv_w2", "std_n_w2", "red_g1_w2", "red_cv_w2", "red_n_w2",
            "ratio_w1", "ratio_w2", "ratio_diff_pct", "ratio", "ratio_window_used", "suppression", "std_g2g1",
            "red_g2g1", "ceiling_gated_pct", "ceil_dna_pct", "twoC_pg_assembly_std", "twoC_pg_literature_std"]
    with open(os.path.join(od, "per_tube.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for row, _, _ in res:
            f.write("\t".join(fmt(row.get(c, "")) for c in cols) + "\n")
    tiers = defaultdict(int)
    for row, _, _ in res:
        tiers[row["tier"]] += 1
    L(f"\nper-tube tiers: {dict(tiers)} (of {len(res)} sample tubes)")

    # ---- 6 summaries
    ok = lambda row: row["tier"] in ("PASS-A", "PASS-B")
    conv = lambda row: row["tier"] not in ("NO-STANDARD",) and not math.isnan(row.get("ratio", np.nan))
    with open(os.path.join(od, "per_stratum.tsv"), "w") as f:
        f.write("date\tdna_channel\tvoltage_set\treference\tref_G1\tn_sample_tubes\tPASS_A\tPASS_B\tFAIL\tn_converged\t"
                "ratio_median_passB\tratio_sd_passB\tratio_cv_pct_passB\tratio_median_converged\tratio_sd_converged\t"
                "ratio_cv_pct_converged\tsuppression_min\tsuppression_max\tstd_cv_w1_median\tred_cv_w1_median\t"
                "twoC_pg_assembly_std_median_passB\n")
        for s in strata:
            rr = [row for row, _, _ in res if (row["date"], row["dna_channel"], row["voltage_set"]) == s]
            if not rr:
                continue
            ref = refs[s]
            pb = [row["ratio"] for row in rr if ok(row)]
            cv_ = [row["ratio"] for row in rr if conv(row)]
            sup = [row["suppression"] for row in rr if conv(row) and not math.isnan(row.get("suppression", np.nan))]
            q = lambda v, fn: fn(v) if len(v) else np.nan
            sd = lambda v: float(np.std(v, ddof=1)) if len(v) > 1 else np.nan
            f.write("\t".join(fmt(x) for x in [
                s[0], s[1], s[2], ref["src"] if ref else "none", ref["g1"] if ref else np.nan, len(rr),
                sum(row["tier"] == "PASS-A" for row in rr), sum(row["tier"] == "PASS-B" for row in rr),
                sum(row["tier"] == "FAIL" for row in rr), len(cv_), q(pb, np.median), sd(pb),
                100 * sd(pb) / np.mean(pb) if len(pb) > 1 else np.nan, q(cv_, np.median), sd(cv_),
                100 * sd(cv_) / np.mean(cv_) if len(cv_) > 1 else np.nan, q(sup, min), q(sup, max),
                q([row["std_cv_w1"] for row in rr if conv(row)], np.nanmedian),
                q([row["red_cv_w1"] for row in rr if conv(row)], np.nanmedian),
                q(pb, np.median) * P["std_2C_pg_assembly"] if len(pb) else np.nan]) + "\n")
    L("\n" + open(os.path.join(od, "per_stratum.tsv")).read())
    # day effect at equal DNA voltage
    byV = defaultdict(dict)
    for row, _, _ in res:
        if conv(row):
            byV[(row["dna_channel"], row["voltage_set"].split("DNA")[1])].setdefault(row["date"], []).append(
                (row["ratio"], ok(row)))
    for k, d in byV.items():
        if len(d) > 1:
            parts = []
            for date, v in sorted(d.items()):
                pv = [x for x, o in v if o]
                parts.append(f"{date}: median {np.median([x for x, _ in v]):.3f} (all converged n {len(v)}), "
                             + (f"{np.median(pv):.3f} (PASS-B+ n {len(pv)})" if pv else "no PASS-B+"))
            L(f"day effect at {k[0]} {k[1]} V: " + " | ".join(parts))
    with open(os.path.join(od, "per_tree.tsv"), "w") as f:
        f.write("date\tvoltage_set\ttree\ttree_class\tn_tubes\tn_chops\tratio_mean_of_chops\tratio_sd_chops\n")
        g = defaultdict(lambda: defaultdict(list))
        for row, _, _ in res:
            if ok(row):
                g[(row["date"], row["voltage_set"], row["tree"], row["tree_class"])][row["chop_id"]].append(row["ratio"])
        for k, ch in sorted(g.items()):
            m = [float(np.mean(v)) for v in ch.values()]
            f.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{k[3]}\t{sum(len(v) for v in ch.values())}\t{len(ch)}\t{np.mean(m):.4f}\t"
                    f"{np.std(m, ddof=1) if len(m) > 1 else float('nan'):.4f}\n")
    rng = np.random.default_rng(int(P["seed"]))
    with open(os.path.join(od, "contrast.tsv"), "w") as f:
        f.write("date\tvoltage_set\tdelta_pct\tci_lo\tci_hi\tn_clone_chops\tn_euploid_trees\tverdict\n")
        for s in strata:
            rr = [row for row, _, _ in res if (row["date"], row["dna_channel"], row["voltage_set"]) == s and ok(row)]
            if not rr:
                continue
            if any(row["tree"] == "UNKNOWN" or row["tree_class"] not in ("clone", "euploid") or row["chop_id"] == "UNKNOWN"
                   for row in rr):
                f.write(f"{s[0]}\t{s[2]}\tNA\tNA\tNA\tNA\tNA\tnot computed: sample sheet has UNKNOWN trees/chops (no contrast)\n")
                continue
            c = fp.contrast(rr, P, rng)
            f.write(f"{s[0]}\t{s[2]}\t" + ("\t".join(fmt(c[k]) for k in ("delta_pct", "ci_lo", "ci_hi", "n_clone_chops",
                                                                         "n_euploid_trees", "verdict")) if c else
                                            "NA\tNA\tNA\tNA\tNA\tnot computed: need 1 clone and >= 2 euploid trees") + "\n")
    L("contrast: " + open(os.path.join(od, "contrast.tsv")).read().replace("\n", " | "))

    # ---- 7 figures
    plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2})
    for s in strata:
        tt = [(row, r, t) for row, r, t in res if (row["date"], row["dna_channel"], row["voltage_set"]) == s]
        stt = [t for t in tubes if t["stratum"] == s and t["sheet"]["role"] == "standard_only"]
        n = len(tt) + len(stt)
        if not n:
            continue
        nc = 4
        nr = math.ceil(n / nc)
        fig, axs = plt.subplots(nr, nc, figsize=(4 * nc, 2.7 * nr), squeeze=False)
        for ax in axs.ravel()[n:]:
            ax.axis("off")
        k = 0
        ref = refs[s]
        for t in stt:
            ax = axs.ravel()[k]; k += 1
            x = t["A"][(t["A"] > 0) & (t["A"] < t["ceil"])]
            ax.hist(x, bins=np.logspace(np.log10(max(x.min(), 10)), np.log10(t["ceil"]), 200), color="#9a9994")
            if ref:
                ax.axvline(ref["g1"], color=INK, lw=0.8, ls="--")
            ax.set_xscale("log"); ax.set_title(f"{t['sheet']['tube']} (standard only)", fontsize=8, loc="left")
        for row, r, t in tt:
            ax = axs.ravel()[k]; k += 1
            if r is None:
                x = t["A"][(t["A"] > 0)]
                ax.hist(np.clip(x, 1, None), bins=np.logspace(0, np.log10(t["ceil"] * 1.05), 200), color="#9a9994")
                ax.axvline(t["ceil"], color="#e34948", lw=0.8)
                ax.set_title(f"{row['tube']}  {row['dna_channel']} {s[2]}\nat ceiling {row['ceil_dna_pct']:.2f}% (no standard)",
                             fontsize=7.5, loc="left")
                ax.set_xscale("log"); continue
            x = r["x"]
            if len(x) < 10:
                ax.set_title(f"{row['tube']}: {len(x)} gated events", fontsize=8, loc="left"); continue
            bins = np.logspace(np.log10(x.min()), np.log10(x.max()), 180)
            ax.hist(x, bins=bins, color="#c3c2b7")
            ax.axvline(P["debris_frac"] * ref["g1"], color=INK2, lw=0.7, ls=":")
            ax.axvspan(P["std_win_lo"] * ref["g1"], P["std_win_hi"] * ref["g1"], color="#2a78d6", alpha=0.08)
            for tag, ls in (("w1", "-"), ("w2", "--")):
                fs_, fr_ = r["fits"].get(tag, (None, None))
                for ff, col in ((fs_, "#2a78d6"), (fr_, "#eb6834")):
                    if ff and ff["converged"]:
                        ax.axvline(ff["mean"], color=col, lw=1 if tag == "w1" else 0.7, ls=ls)
            ax.set_xscale("log")
            ax.set_title(f"{row['tube']}  [{row['tier']}]\nratio {fmt(row['ratio'])[:5]}  n {row['std_n_w1']}/{row['red_n_w1']}  "
                         f"CV {row['std_cv_w1']:.1f}/{row['red_cv_w1']:.1f}%", fontsize=7.5, loc="left")
        fig.suptitle(f"{a.run}  {s[0]}  {s[1]}  {s[2]}  (blue: standard window & G1; orange: redwood G1; solid = +-"
                     f"{int(P['fit_w1']*100)}% fit, dashed = +-{int(P['fit_w2']*100)}%; dotted = debris cut)", fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(fd, f"hist_{s[0].replace('-','')}_{s[2].replace('/','_')}.png"), dpi=130)
        plt.close(fig)
    rr = [row for row, _, _ in res if conv(row)]
    if rr:
        lab = sorted({(row["date"], row["voltage_set"]) for row in rr})
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw=dict(width_ratios=[1.4, 1]))
        jit = np.random.default_rng(1)
        for row in rr:
            xi = lab.index((row["date"], row["voltage_set"])) + jit.uniform(-0.15, 0.15)
            ax1.scatter(xi, row["ratio"], s=34, marker=CLASSMK.get(row["tree_class"], "D"),
                        facecolor=TIERCOL[row["tier"]] if row["tier"] != "FAIL" else "none",
                        edgecolor=TIERCOL[row["tier"]], lw=1.2, zorder=3)
            if not math.isnan(row.get("suppression", np.nan)):
                ax2.scatter(row["suppression"], row["ratio"], s=30, marker=CLASSMK.get(row["tree_class"], "D"),
                            facecolor=TIERCOL[row["tier"]] if row["tier"] != "FAIL" else "none",
                            edgecolor=TIERCOL[row["tier"]], lw=1.2)
        ax1.set_xticks(range(len(lab))); ax1.set_xticklabels([f"{d}\n{v}" for d, v in lab], fontsize=7.5)
        ax1.set_ylabel("redwood G1 / in-tube standard G1"); ax1.yaxis.grid(True, color=GRID, lw=0.5)
        ax1.set_title("ratio by day and voltage set (filled = PASS-B or better)", fontsize=9, loc="left")
        ax2.set_xlabel("suppression (in-tube / standard-only G1)"); ax2.set_ylabel("ratio")
        ax2.set_title("suppression vs ratio", fontsize=9, loc="left"); ax2.yaxis.grid(True, color=GRID, lw=0.5)
        h = [plt.Line2D([], [], ls="", marker="o", color=c, label=k) for k, c in TIERCOL.items()]
        h += [plt.Line2D([], [], ls="", marker=m, color=INK2, mfc="none", label=f"tree: {k}") for k, m in CLASSMK.items()
              if k != "standard"]
        ax1.legend(handles=h, frameon=False, fontsize=7, ncol=2)
        for ax in (ax1, ax2):
            ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout(); fig.savefig(os.path.join(fd, "ratio_by_stratum_and_suppression.png"), dpi=160); plt.close(fig)
    L(f"\noutputs: {od}; figures: {fd}")


if __name__ == "__main__":
    main()
