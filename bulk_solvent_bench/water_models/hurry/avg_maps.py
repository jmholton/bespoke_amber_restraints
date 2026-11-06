# frame-weighted average of per-iteration solvent maps (bfac_iter.csh: <dir>/avg.map, <dir>/nframes) -> FPART PHIPART
# usage (phenix.python, octamus1): avg_maps.py out.mtz dir1 dir2 ... [reso=0.965]
import sys, os, subprocess
import numpy as np
import gemmi
GEMMI = "/programs/ccp4-9/bin/gemmi"
out = sys.argv[1]
dirs = [a for a in sys.argv[2:] if "=" not in a]
opts = dict(a.split("=") for a in sys.argv[2:] if "=" in a)
reso = opts.get("reso", "0.965")
total, nf, m = None, 0, None
for d in dirs:
    n = int(open(os.path.join(d, "nframes")).read())
    m = gemmi.read_ccp4_map(os.path.join(d, "avg.map"))
    g = np.array(m.grid, copy=True) * n
    total = g if total is None else total + g
    nf += n
m.grid.array[:] = (total / nf).astype(np.float32)
m.update_ccp4_header()
avgmap = out[:-4] + "_avg.map"
m.write_ccp4_map(avgmap)
subprocess.run([GEMMI, "map2sf", "--dmin=" + reso, avgmap, out + ".tmp.mtz", "FPART", "PHIPART"], check=True,
               stdout=subprocess.DEVNULL)
os.replace(out + ".tmp.mtz", out)
print("%s: %d iterations, %d frames, mean density %.5f e/A^3" % (out, len(dirs), nf, total.mean() / nf))
