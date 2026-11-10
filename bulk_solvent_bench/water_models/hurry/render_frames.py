# average solvent density of a set of PDB frames, rendered as nc2mtz_gemmi.com renders each frame: every atom at
# occupancy 1 and one uniform B, gemmi sfcalc in the crystal's space group (which lays the MD cell's symmetry copies
# on top of each other), summed, scaled by 1/(nsym x nframes), gemmi map2sf -> FPART PHIPART (electron units).
# usage (phenix.python, octamus1): render_frames.py out.mtz frame_list.txt [B=10] [reso=0.965] [nproc=48]
import sys, os, subprocess, tempfile
from multiprocessing import Pool
import numpy as np
import gemmi
CELL = "CRYST1   45.900   40.700   30.100  90.00  90.00  90.00 P 21 21 21"
NSYM = 4
GEMMI = "/programs/ccp4-9/bin/gemmi"
out, listfile = sys.argv[1], sys.argv[2]
opts = dict(a.split("=") for a in sys.argv[3:])
B = float(opts.get("B", 10))
reso = opts.get("reso", "0.965")
nproc = int(opts.get("nproc", 48))
frames = [l.strip() for l in open(listfile) if l.strip()]
tmp = tempfile.mkdtemp(prefix="render_", dir="/dev/shm")

def render(i_f):
    i, f = i_f
    pdb = os.path.join(tmp, "%d.pdb" % i)
    with open(pdb, "w") as w:
        w.write(CELL + "\n")
        n = 0
        for line in open(f):
            if not line.startswith(("ATOM", "HETATM")):
                continue
            n += 1
            el = line[76:78].strip() or line[12:14].strip()
            # the line nc2mtz's per-frame job writes: element in the name field, occ, B, element in cols 77-78
            w.write("ATOM  %5d %-2s%40s%6.2f%6.2f%12s%12d\n" % (n % 100000, el, line[14:54], 1.0, B, el, n))
        w.write("END\n")
    mp, mtz = pdb[:-4] + ".map", pdb[:-4] + ".mtz"
    subprocess.run([GEMMI, "sfcalc", "--dmin=" + reso, "--rate=2", "--write-map=" + mp, "--to-mtz=" + mtz, pdb],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    g = np.array(gemmi.read_ccp4_map(mp).grid, copy=True)
    for x in (pdb, mp, mtz):
        os.remove(x)
    return g

total, n = None, 0
with Pool(nproc) as pool:
    for g in pool.imap_unordered(render, list(enumerate(frames)), chunksize=4):
        total = g if total is None else total + g
        n += 1
ref = os.path.join(tmp, "ref.pdb")
mp0 = os.path.join(tmp, "grid.map")
with open(ref, "w") as w:
    w.write(CELL + "\nATOM      1 O    HOH A   1       1.000   1.000   1.000  1.00 10.00           O\nEND\n")
subprocess.run([GEMMI, "sfcalc", "--dmin=" + reso, "--rate=2", "--write-map=" + mp0, "--to-mtz=" + mp0 + ".mtz", ref],
               check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
m = gemmi.read_ccp4_map(mp0)
assert m.grid.shape == total.shape
m.grid.array[:] = (total / (NSYM * n)).astype(np.float32)
m.update_ccp4_header()
avgmap = out[:-4] + "_avg.map"
m.write_ccp4_map(avgmap)
subprocess.run([GEMMI, "map2sf", "--dmin=" + reso, avgmap, out + ".tmp.mtz", "FPART", "PHIPART"], check=True,
               stdout=subprocess.DEVNULL)
os.replace(out + ".tmp.mtz", out)
for x in os.listdir(tmp):
    os.remove(os.path.join(tmp, x))
os.rmdir(tmp)
print("%s: %d frames, B %g, mean density %.5f e/A^3" % (out, n, B, total.mean() / (NSYM * n)))
