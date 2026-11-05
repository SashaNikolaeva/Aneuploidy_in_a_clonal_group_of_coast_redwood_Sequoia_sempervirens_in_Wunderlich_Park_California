# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Strict-gate genome sizes of the same tubes (gate comparison behind the Discussion sentence
# on variation between analyses). Not the script behind the plotted values.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
import fcsparser, numpy as np, glob, os, csv
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = "/Users/sashanikolaeva/Desktop/Flow_data/Experiment_002"
OUT = "/private/tmp/claude-501/-Users-sashanikolaeva/b0e7c8c5-0b79-4ece-a7f5-e908d07d0ee1/scratchpad"
files = sorted(glob.glob(os.path.join(DIR, "*.fcs")))
A, H = "BUV 496-A", "BUV 496-H"
# The G1 peak = 2C = total DNA content of the resting somatic nucleus (makes NO
# assumption about the sample's ploidy). SEGI is a well-characterized diploid standard:
# its G1 nucleus holds 2C = 2 x the 8.125 Gb haploid assembly = 16.25 Gb.
SEGI_2C_GBP = 16.25
GBP_PER_PG  = 0.978
CVMAX = 22.0

def label(fn): return os.path.basename(fn).replace("Specimen_001_","").replace(".fcs","")
def gauss(x,a,mu,sig): return a*np.exp(-0.5*((x-mu)/sig)**2)

def voltage(fn):
    meta=fcsparser.parse(fn, meta_data_only=True); n=int(meta['$PAR'])
    for i in range(1,n+1):
        if meta.get(f'$P{i}N')==A: return meta.get(f'$P{i}V')
    return "?"

def load_gated(fn, debris=10000, lo=0.7, hi=1.4):
    meta, df = fcsparser.parse(fn, reformat_meta=True)
    a = df[A].values.astype(float); h = df[H].values.astype(float)
    g = np.isfinite(a)&np.isfinite(h)&(a>0)&(h>0); a,h=a[g],h[g]
    sel=a>debris; a2,h2=a[sel],h[sel]
    if len(a2)>50:
        r=a2/h2; m=np.median(r); a2=a2[(r>m*lo)&(r<m*hi)]   # singlet band (doublet gate)
    satpct=100.0*np.mean(a>=255000)
    return a2, satpct

def fit_peak(a, center, half=9000):
    seg=a[(a>=center-half)&(a<=center+half)]
    if len(seg)<20: return None
    cnt,edg=np.histogram(seg,bins=40,range=(center-half,center+half)); ctr=0.5*(edg[:-1]+edg[1:])
    try:
        popt,_=curve_fit(gauss,ctr,cnt,p0=[cnt.max(),center,half/3],maxfev=5000)
        mu,sig=popt[1],abs(popt[2])
        if not (center-half<mu<center+half): mu,sig=np.median(seg),seg.std()
    except Exception:
        mu,sig=np.median(seg),seg.std()
    return mu, sig/mu*100.0, len(seg)

def two_peaks(a):
    top=min(a.max(),262143)
    cnt,edg=np.histogram(a,bins=400,range=(10000,top*1.02)); ctr=0.5*(edg[:-1]+edg[1:])
    sc=gaussian_filter1d(cnt.astype(float),2)
    pk,props=find_peaks(sc,prominence=max(sc.max()*0.06,3),distance=6)
    cand=[(ctr[p],props['prominences'][k]) for k,p in enumerate(pk) if ctr[p]<255000]
    fits=[]
    for pos,_ in sorted(cand,key=lambda t:-t[1])[:4]:   # fit up to 4 most prominent
        f=fit_peak(a,pos)
        if f and f[1]<=CVMAX: fits.append(f)
    fits=sorted(fits,key=lambda t:t[0])                 # by position
    return fits

def classify(name):
    base=name.split("_")[0]
    if base=="SEGI": return "SEGI control (standard only)"
    has_C = any(c in base for c in ["C","c"]) and base!="SEGI"
    has_NR= "NR" in base
    if has_NR and ("c2" in base or "C" in base.replace("NR","")):
        return "normal redwood + aneuploid clone"
    if has_NR: return "normal redwood (not clonal)"
    if base.startswith(("R","Root","Roof")): return "root (clonal, aneuploid redwood)"
    return "aneuploid redwood clone"

def assign(a, stype):
    peaks=two_peaks(a) if len(a)>=60 else []
    segi=sample=None; note=""
    if "control" in stype:
        if peaks: segi=peaks[-1]; note="reference position only"
    elif len(peaks)>=2:
        segi=peaks[0]; sample=peaks[-1]      # lower=SEGI standard, upper=redwood/clone
    elif len(peaks)==1:
        note="only one peak resolved"
    else:
        note="no clean nuclei peaks"
    return peaks, segi, sample, note

rows=[]; fig,axes=plt.subplots(6,3,figsize=(16,22)); axes=axes.ravel()
for ax,fn in zip(axes,files):
    name=label(fn); V=voltage(fn); stype=classify(name)
    a,satpct=load_gated(fn)                        # strict gate (0.7-1.4x, debris>10000)
    peaks,segi,sample,note=assign(a,stype)
    gate="strict"
    # Adaptive rescue: if a non-control tube can't resolve BOTH peaks, retry a
    # RELAXED gate (0.6-1.7x, debris>5000). The doublet A/B test showed this
    # biases 2C only ~1% high, so it is a safe last resort -- but flag it.
    if ("control" not in stype) and not (segi and sample):
        a2,sat2=load_gated(fn,debris=5000,lo=0.6,hi=1.7)
        p2,s2,m2,n2=assign(a2,stype)
        if s2 and m2:
            a,satpct,peaks,segi,sample=a2,sat2,p2,s2,m2
            gate="relaxed"; note="RELAXED gate (2C may read ~1% high)"

    top=min(a.max(),262143) if len(a) else 1
    cnt,edg=np.histogram(a,bins=300,range=(10000,top*1.02)); ctr=0.5*(edg[:-1]+edg[1:])
    sc=gaussian_filter1d(cnt.astype(float),2); ax.plot(ctr,sc,lw=0.8)
    ax.set_title(f"{name}  {V}V (n={len(a)}, {gate})",fontsize=9); ax.set_xlabel(A)

    for i,(mu,cv,nn) in enumerate(peaks):
        role = "SEGI" if (segi and mu==segi[0]) else ("redwood/clone" if (sample and mu==sample[0]) else "")
        col = 'tab:red' if role=="SEGI" else ('tab:purple' if role=="redwood/clone" else 'gray')
        ax.axvline(mu,color=col,ls=':',lw=0.9)
        ax.text(mu,sc.max()*0.88,f"{role}\n{mu:.0f}\nCV{cv:.0f}%",fontsize=6,ha='center')

    ratio=gb2c=pg2c=np.nan
    if segi and sample:
        ratio=sample[0]/segi[0]            # ratio of G1 (2C) peaks = ratio of total DNA content
        gb2c=ratio*SEGI_2C_GBP; pg2c=gb2c/GBP_PER_PG
        if not (2.3<=ratio<=4.2): note=(note+" | ratio atypical").strip(" |")
    if satpct>0.5: note=(note+f" | {satpct:.1f}% near scale ceiling").strip(" |")

    cvs=[v[1] for v in [segi,sample] if v]
    qual = "--" if not cvs else ("excellent" if max(cvs)<8 else "good" if max(cvs)<12 else "fair" if max(cvs)<18 else "poor")
    rows.append(dict(sample=name, voltage=V, sample_type=stype, gate=gate, quality=qual, n_nuclei=len(a),
        SEGI_mean=round(segi[0]) if segi else "", SEGI_CV=round(segi[1],1) if segi else "",
        sample_mean=round(sample[0]) if sample else "", sample_CV=round(sample[1],1) if sample else "",
        ratio=round(ratio,3) if ratio==ratio else "",
        DNA_2C_Gb=round(gb2c,2) if gb2c==gb2c else "",
        DNA_2C_pg=round(pg2c,2) if pg2c==pg2c else "", note=note))

for ax in axes[len(files):]: ax.axis("off")
plt.tight_layout()
cols=["sample","sample_type","gate","voltage","quality","n_nuclei","SEGI_mean","SEGI_CV","sample_mean","sample_CV",
      "ratio","DNA_2C_Gb","DNA_2C_pg","note"]
for d in (OUT,DIR):
    if not os.path.isdir(d): continue
    with open(os.path.join(d,"genome_sizes.csv"),"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=cols); w.writeheader(); w.writerows(rows)
    plt.savefig(os.path.join(d,"genome_sizes.png"),dpi=95)

h=f"{'sample':24s}{'type':34s}{'V':>4s}{'SEGI':>7s}{'smpl':>8s}{'ratio':>6s}{'2C_Gb':>7s}{'2C_pg':>7s}  qual"
print(h); print("-"*len(h))
for r in rows:
    print(f"{r['sample']:24s}{r['sample_type']:34s}{str(r['voltage']):>4s}{str(r['SEGI_mean']):>7s}"
          f"{str(r['sample_mean']):>8s}{str(r['ratio']):>6s}{str(r['DNA_2C_Gb']):>7s}{str(r['DNA_2C_pg']):>7s}  {r['quality']}  {r['note']}")
print("\nwrote genome_sizes.csv + genome_sizes.png to data folder")
