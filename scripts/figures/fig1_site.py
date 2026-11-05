#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Figure 1 (range map, site map, historical and present imagery, site photos).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""New Figure 1: where the Wunderlich Park clone is, how old it is at least, and what it looks like.

Panels
  A  range map: California/Oregon state lines, the survey trees (canonical 224), the clone (star). (The Little 1971
     range polygon was removed 2026-10-26 at Sasha's request.)
  B  site map: the 9 ramets, WPR and C1 at their V5 GPS points (Data/Coordinates/Updated_coordinates.V5.csv via
     the tree inventory), local metres east/north of the ramet centroid, scale bar, north arrow. GPS points that
     coincide (13_WP1 = 22_WPd2, 14_WP3 = 31_WP2_1) are drawn as separate points SPREAD m apart and marked *.
  C  historical imagery, Google Earth 26 Sep 1948 (earliest layer; Sasha's screenshot, clone outlined by her)
  D  present-day imagery, Google Earth 5 Mar 2026 (same view); both cropped to CROP_GE to remove the map controls
  E, F  site photos (Sasha, Results/clone/Photos)                            <- PHOTOS below
Output: clone_paper_figures/site/fig1_site.{png,pdf}
Run: ~/miniconda3/envs/quast/bin/python3 -I scripts/figures/fig1_site.py
"""
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402
from matplotlib.patches import Polygon as MPoly  # noqa: E402
from PIL import Image, ImageDraw, ImageOps  # noqa: E402

S1 = "/space/s1/sashanikolaeva"
POPMAP = f"{S1}/Scripts/GitHub/redwood-reference-genome/scripts/popmap/build/src"
INV = f"{S1}/Scripts/GitHub/redwood-reference-genome/results/popmap/redwood_tree_inventory.tsv"
PH = f"{S1}/Results/clone/Photos"
PHOTOS = [f"{PH}/31049F29-E1E0-45B8-9AAD-ACE1287AC353.JPG",   # canopy opening at the centre of the ring
          f"{PH}/F710A3F3-51AA-4D76-9EAC-0AF4A4599374.JPG"]   # leaning ramets and fallen stems
# Sasha's screenshots (2872 x 1156; macOS names contain a narrow no-break space, hence glob)
HIST = glob.glob(f"{PH}/Screenshot 2026-10-26 at 3.27.10*PM.png")[0]   # date stamp 9/26/1948
NOW = glob.glob(f"{PH}/Screenshot 2026-10-26 at 3.30.09*PM.png")[0]    # date stamp 3/5/2026
CROP_GE = (970, 50, 1970, 1050)   # left, top, right, bottom: clear of the date pill, minimap and buttons
# Sasha's 2026 outline has a gap at its NW corner; close it with a segment of the same colour and width between the
# two open ends (screenshot pixel coordinates, measured from the yellow mask; Sasha asked 2026-10-26)
NOW_GAP = [(1362, 425), (1399, 397)]
OUTLINE_RGB, OUTLINE_W = (238, 190, 75), 5
# north is up in both screenshots (Google Earth compass needle points up). Metres per screenshot pixel for a scale bar:
# None until Sasha measures a distance with the Google Earth ruler (not derivable from the screenshot itself)
GE_M_PER_PX = {"C": None, "D": None}
SCALE_M = 50
HIST_LABEL, NOW_LABEL = "26 September 1948 (Google Earth historical imagery)", "5 March 2026 (Google Earth)"
OUT = f"{S1}/Scripts/GitHub/clone_paper_figures/site"
os.makedirs(OUT, exist_ok=True)

CLONE = ["13_WP1", "6_WP5", "31_c3", "31_WP2_1", "28_WP7", "28_WP6", "14_WP3", "22_WPd2", "33_c6"]
GREEN, ORANGE, LORANGE = "#1b7837", "#e08214", "#fdb863"


def load_js(path, var):
    s = open(path).read()
    return json.loads(s[s.index("{"):s.rindex("}") + 1])


def topo_polys(topo, obj):
    """Decode a TopoJSON GeometryCollection into a list of rings (lon, lat arrays)."""
    sx, sy = topo["transform"]["scale"]
    tx, ty = topo["transform"]["translate"]
    arcs = []
    for a in topo["arcs"]:
        q = np.cumsum(np.array(a, dtype=float), axis=0)
        arcs.append(np.column_stack([q[:, 0] * sx + tx, q[:, 1] * sy + ty]))

    def ring(ix):
        pts = []
        for i in ix:
            a = arcs[i] if i >= 0 else arcs[~i][::-1]
            pts.append(a if not pts else a[1:])
        return np.vstack(pts)

    out = []
    for g in topo["objects"][obj]["geometries"]:
        polys = g["arcs"] if g["type"] == "MultiPolygon" else [g["arcs"]] if g["type"] == "Polygon" else []
        for p in polys:
            out.extend(ring(r) for r in p)
    return out


inv = pd.read_csv(INV, sep="\t", low_memory=False).set_index("sample_id")
survey = inv[inv.in_canonical_224 == True]  # noqa: E712
c_lat, c_lon = inv.loc[CLONE, "lat"].mean(), inv.loc[CLONE, "lon"].mean()

fig = plt.figure(figsize=(15, 11.2))
gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.78], wspace=0.18, hspace=0.16)

# ---- A: range map
ax = fig.add_subplot(gs[0, 0])
topo = load_js(f"{POPMAP}/topo.js", "TOPO")
k = np.cos(np.radians(39.5))
for r in topo_polys(topo, "states"):
    if r[:, 0].max() < -125.5 or r[:, 0].min() > -116 or r[:, 1].max() < 35 or r[:, 1].min() > 43.5:
        continue
    ax.add_patch(MPoly(np.column_stack([r[:, 0] * k, r[:, 1]]), closed=True, fc="#f4f1ea", ec="#9e9a90", lw=0.6))
ax.scatter(survey.lon * k, survey.lat, s=5, c="0.25", lw=0, zorder=3, label="survey trees (n = 224)")
ax.scatter([c_lon * k], [c_lat], marker="*", s=260, c=GREEN, ec="k", lw=0.7, zorder=4, label="Wunderlich Park clone")
ax.set_xlim(-125.0 * k, -119.6 * k)
ax.set_ylim(35.4, 42.6)
ax.set_aspect("equal")
ax.set_facecolor("#dbe9f2")
ax.set_xticks([])
ax.set_yticks([])
# legend in the empty interior (upper right), narrow two-line labels so it covers no trees
# (redwood range polygon removed at Sasha's request, 2026-10-26)
hand = [plt.Line2D([], [], ls="", marker="o", ms=3, color="0.25", label="survey trees\n(n = 224)"),
        plt.Line2D([], [], ls="", marker="*", ms=12, mfc=GREEN, mec="k", mew=0.6, label="Wunderlich\nPark clone")]
ax.legend(handles=hand, loc="upper right", fontsize=7.5, frameon=True, framealpha=0.95, borderpad=0.6, handletextpad=0.4,
          labelspacing=0.9)
ax.text(-120.9 * k, 36.6, "California", fontsize=9, color="0.45", style="italic")
ax.annotate("San Francisco\nBay", (-122.35 * k, 37.7), xytext=(-121.3 * k, 38.0), fontsize=7.5, color="0.35",
            arrowprops=dict(arrowstyle="-", color="0.5", lw=0.6))
ax.set_title("A", loc="left", fontweight="bold", fontsize=15)

# ---- B: site map (metres)
ax = fig.add_subplot(gs[0, 1])
pts = CLONE + ["22_WPR", "16_C1"]
lat, lon = inv.loc[pts, "lat"], inv.loc[pts, "lon"]
x = (lon - c_lon) * 111320 * np.cos(np.radians(c_lat))
y = (lat - c_lat) * 110574
xy = pd.DataFrame(dict(x=x, y=y))
# stems sharing a GPS fix (13_WP1 = 22_WPd2, 31_WP2_1 = 14_WP3) are drawn as separate points, spread east-west by
# SPREAD m around the shared fix (Sasha, 2026-10-26); the true positions are identical
SPREAD = 7.0
xy_true = xy.copy()
done = set()
for s in CLONE:
    if s in done:
        continue
    same = [t for t in CLONE if abs(xy_true.x[t] - xy_true.x[s]) < 0.5 and abs(xy_true.y[t] - xy_true.y[s]) < 0.5]
    done.update(same)
    if len(same) > 1:
        for i, t in enumerate(same):
            xy.loc[t, "x"] = xy_true.x[t] + (i - (len(same) - 1) / 2) * SPREAD
shared = {t for t in CLONE if (abs(xy_true.x - xy_true.x[t]) < 0.5).sum() > 1 and
          ((abs(xy_true.x - xy_true.x[t]) < 0.5) & (abs(xy_true.y - xy_true.y[t]) < 0.5)).sum() > 1}
cl = xy.loc[CLONE]
ax.scatter(cl.x, cl.y, s=110, c=GREEN, ec="k", lw=0.6, zorder=3, label="clone ramets (n = 9)")
for s in CLONE:
    lab = s.split("_", 1)[1] + (" *" if s in shared else "")
    left = s in shared and xy.x[s] < xy_true.x[s]
    ax.annotate(lab, (xy.x[s], xy.y[s]), xytext=(-7 if left else 7, 6), textcoords="offset points", fontsize=9,
                color=GREEN, ha="right" if left else "left")
s = "22_WPR"
ax.scatter(xy.x[s], xy.y[s], s=95, c=ORANGE, ec="k", lw=0.6, zorder=4, marker="D")
ax.annotate("WPR", (xy.x[s], xy.y[s]), xytext=(8, -12), textcoords="offset points", fontsize=9.5)
pad = 22
ax.set_aspect("equal")
ax.set_xlim(cl.x.min() - pad - 15, cl.x.max() + pad)
ax.set_ylim(cl.y.min() - pad - 10, cl.y.max() + pad + 12)
# C1 lies outside the panel: arrow toward it
dx, dy = xy.x["16_C1"], xy.y["16_C1"]
# label distance: haversine from the clone centroid (236 m, as in Methods); the local projection gives 235 m
r_ = np.pi / 180
la1, lo1, la2, lo2 = c_lat * r_, c_lon * r_, inv.loc["16_C1", "lat"] * r_, inv.loc["16_C1", "lon"] * r_
dist = float(2 * 6371000 * np.arcsin(np.sqrt(np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2)))
ux, uy = dx / np.hypot(dx, dy), dy / np.hypot(dx, dy)
x1, y1 = cl.x.min() - pad - 10, cl.y.max() + pad + 2
ax.annotate("", xy=(x1, y1), xytext=(x1 - 30 * ux, y1 - 30 * uy),
            arrowprops=dict(arrowstyle="-|>", color=LORANGE, lw=2))
ax.text(x1 + 6, y1 + 1, f"C1: {dist:.0f} m to the NW", fontsize=9, va="top", color="#b35806")
x0, y0 = cl.x.min() - pad - 10, cl.y.min() - pad - 2
ax.plot([x0, x0 + 25], [y0, y0], color="k", lw=2.5)
ax.text(x0 + 12.5, y0 + 3, "25 m", ha="center", fontsize=8.5)
xn, yn = cl.x.max() + pad - 8, cl.y.min() - pad + 2
ax.annotate("N", (xn, yn + 20), xytext=(xn, yn), ha="center", fontsize=10, fontweight="bold",
            arrowprops=dict(arrowstyle="-|>", color="k", lw=1.2))
ax.set_xlabel("metres east of the clone centre")
ax.set_ylabel("metres north of the clone centre")
ax.text(0.55, 0.02, f"* same GPS fix, drawn {SPREAD:.0f} m apart\n(handheld GPS under canopy, ~5-10 m)", transform=ax.transAxes, va="bottom", ha="center", fontsize=7.5,
        color="0.4")
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("B", loc="left", fontweight="bold", fontsize=15)


def image_panel(spec, path, letter, label, crop=None, close_gap=None, north=False):
    ax = fig.add_subplot(spec)
    if os.path.exists(path):
        im = ImageOps.exif_transpose(Image.open(path)).convert("RGB")
        if close_gap:
            ImageDraw.Draw(im).line(close_gap, fill=OUTLINE_RGB, width=OUTLINE_W, joint="curve")
        if crop:
            im = im.crop(crop)
        w0 = im.width
        im.thumbnail((1800, 1800))
        ax.imshow(im)
        if north:  # north arrow, white with a dark edge so it reads on the imagery
            ax.annotate("N", xy=(0.93, 0.95), xytext=(0.93, 0.80), xycoords="axes fraction", ha="center", va="top",
                        fontsize=11, fontweight="bold", color="white",
                        arrowprops=dict(arrowstyle="-|>", color="white", lw=2, mutation_scale=16),
                        path_effects=[pe.withStroke(linewidth=2.5, foreground="black")])
        mpp = GE_M_PER_PX.get(letter)
        if mpp:  # scale bar of SCALE_M metres, in displayed pixels
            L = SCALE_M / mpp * im.width / w0
            x0, y0 = 0.05 * im.width, 0.94 * im.height
            ax.plot([x0, x0 + L], [y0, y0], color="white", lw=4, solid_capstyle="butt",
                    path_effects=[pe.withStroke(linewidth=6, foreground="black")])
            ax.text(x0 + L / 2, y0 - 0.025 * im.height, f"{SCALE_M} m", ha="center", va="bottom", fontsize=9,
                    color="white", path_effects=[pe.withStroke(linewidth=2.5, foreground="black")])
    else:
        ax.set_facecolor("#eeeeee")
        ax.text(0.5, 0.5, f"{label}\n(file not yet uploaded)", ha="center", va="center", fontsize=10, color="0.4",
                transform=ax.transAxes)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(letter, loc="left", fontweight="bold", fontsize=15)
    if label:
        ax.set_xlabel(label, fontsize=9)


image_panel(gs[1, 0], HIST, "C", HIST_LABEL, CROP_GE, north=True)
image_panel(gs[1, 1], NOW, "D", NOW_LABEL, CROP_GE, close_gap=NOW_GAP, north=True)
image_panel(gs[0, 2], PHOTOS[0], "E", "canopy opening at the centre of the ring")
image_panel(gs[1, 2], PHOTOS[1], "F", "leaning ramets and fallen stems")

fig.savefig(f"{OUT}/fig1_site.png", dpi=200, bbox_inches="tight")
fig.savefig(f"{OUT}/fig1_site.pdf", bbox_inches="tight", dpi=300)
print("clone centre", round(c_lat, 5), round(c_lon, 5))
print(xy.round(1).to_string())
