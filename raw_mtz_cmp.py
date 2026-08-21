#!/usr/bin/env ccp4-python
"""Compare two sfcalc_gpu_collapse MTZ files by reading their data blocks directly.

The C++ exe's MTZ header is rejected by gemmi and CCP4 mtzdmp (its SYMINF record
is missing the point-group token), so nothing can open these files normally.  The
data block itself is unambiguous: "MTZ " + int header_loc (1-indexed word count)
+ 4-byte machine stamp, then nrefl*ncol float32 in row-major order.

usage:  raw_mtz_cmp.py a.mtz b.mtz [ncol]
          ncol 5 (default) = H K L FC PHIC   (outmtz)
          ncol 4           = H K L I         (outI)
"""
import sys, struct, numpy as np

def read(path, ncol):
    b = open(path, 'rb').read()
    if b[:4] != b'MTZ ':
        sys.exit("%s is not an MTZ file" % path)
    nwords = struct.unpack('<i', b[4:8])[0] - 4
    if nwords % ncol:
        sys.exit("%s: %d data words is not a multiple of ncol=%d" % (path, nwords, ncol))
    a = np.frombuffer(b, dtype='<f4', count=nwords, offset=12).reshape(-1, ncol)
    return {tuple(r): a[i, 3:] for i, r in enumerate(a[:, :3].astype(int))}

def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    ncol = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    d1, d2 = read(sys.argv[1], ncol), read(sys.argv[2], ncol)
    keys = sorted(set(d1) & set(d2))
    if not keys:
        sys.exit("no reflections in common")
    F1 = np.array([d1[k][0] for k in keys])
    F2 = np.array([d2[k][0] for k in keys])
    # a global scale is absorbed downstream, so scale before judging
    k = (F1 * F2).sum() / (F2 * F2).sum()
    out = "n=%d  scale=%.5f  R=%.5f  CC=%.6f" % (
        len(keys), k, np.abs(F1 - k * F2).sum() / np.abs(F1).sum(),
        np.corrcoef(F1, F2)[0, 1])
    if ncol >= 5:   # phased: also compare as complex numbers
        C1 = F1 * np.exp(1j * np.deg2rad([d1[key][1] for key in keys]))
        C2 = F2 * np.exp(1j * np.deg2rad([d2[key][1] for key in keys]))
        out += "  R_complex=%.5f" % (np.abs(C1 - k * C2).sum() / np.abs(C1).sum())
    print("%s vs %s: %s" % (sys.argv[1], sys.argv[2], out))

main()
