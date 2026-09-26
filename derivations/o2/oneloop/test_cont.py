"""DIS kinematics: MCFM with the standard -i0 continuation vs with the
imaginary parts of lnrat dropped, both fitted to DISENT's photon
non-factorising part (2 LEIV - EMSQ/2 ERTV)/EMSQ plus tree terms."""
import sys, math, numpy as np, importlib
sys.path.insert(0, "oneloop")
import bdk
exec(open("oneloop/fit_photon.py").read().split("def fit(")[0].replace("from bdk import spinoru, Amp", "from bdk import spinoru, Amp"))
def fitres(rows):
    out = []
    for tgt in ("I51", "I52"):
        X = np.array([[r["NX"], r["NY"], r["T"]] + [r["T"]*l for l in r["L"]] +
                      [r["T"]*r["L"][i]*r["L"][j] for i in range(3) for j in range(i, 3)] for r in rows])
        y = np.array([r[tgt] for r in rows])
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        out.append((tgt, np.round(c[:2], 8), np.max(np.abs(X @ c - y)/np.abs(y))))
    return out
print("standard continuation:", fitres(rows))
