# how flat each solvent map is away from the protein: mean and sd of the density by distance from the nearest protein
# heavy atom of the example model (symmetry included), and complex CC with the 2019 map by resolution shell.
# usage: phenix.python bulk_flatness.py model.pdb name=map.mtz ...
import sys
import numpy as np
import gemmi
st = gemmi.read_structure(sys.argv[1])
st.setup_entities()
maps = [a.split("=", 1) for a in sys.argv[2:]]
def fmap(f):
    m = gemmi.read_mtz_file(f)
    return m
ref = fmap(maps[0][1])
grid0 = ref.transform_f_phi_to_map("FPART", "PHIPART", sample_rate=3)
# distance to the nearest protein heavy atom on that grid, in 0.5 A steps up to 8 A
dist = np.full(np.array(grid0, copy=False).shape, 99.0, dtype=np.float32)
for r in np.arange(8.0, 0.49, -0.5):
    mk = gemmi.FloatGrid(*dist.shape)
    mk.set_unit_cell(grid0.unit_cell)
    mk.spacegroup = grid0.spacegroup
    for ch in st[0]:
        for res in ch:
            if res.name == "HOH":
                continue
            for a in res:
                if a.element.name not in ("H", "D"):
                    mk.set_points_around(a.pos, float(r), 1.0)
    mk.symmetrize_max()
    dist[np.array(mk, copy=False) > 0] = r
bins = [(2.5, 3.5), (3.5, 5.0), (5.0, 6.5), (6.5, 99)]
print("fraction of cell:", " ".join("%.1f-%s A %.3f" % (lo, hi if hi < 99 else "", ((dist > lo) & (dist <= hi)).mean()) for lo, hi in bins))
md = None
for name, f in maps:
    m = fmap(f)
    g = np.array(m.transform_f_phi_to_map("FPART", "PHIPART", exact_size=list(dist.shape)), copy=False)
    row = []
    for lo, hi in bins:
        s = (dist > lo) & (dist <= hi)
        row.append("%.3f/%.3f" % (g[s].mean(), g[s].std()))
    # sd/mean in the farthest bin: how spiky the bulk is
    s = dist > 5.0
    print("%-12s mean/sd by distance: %s   bulk sd/mean (>5 A) %.2f" % (name, "  ".join(row), g[s].std() / g[s].mean()))
