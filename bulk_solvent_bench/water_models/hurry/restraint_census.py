# usage: python3 restraint_census.py rundir first last step   (on octamus1; rundir in 1aho_hurry_test)
# Which atoms a bespoke optimize_weights run restrained during its MD, and how hard, averaged over iterations
# first..last.  The MD input amber_<itr>.in carries the positional restraints (ntr=1) as groups: "Specific <atom>
# restraints on N atoms in range: lo-hi", the group's weight (kcal/mol/A^2) on the next line, FIND <atom name>, then
# SEARCH RES ranges.  Residue names come from the run's padded.parm7: the MD topology is padded.parm7 cut at the end,
# so residue numbers agree.  Water and ion totals are the residues present in that iteration (nwaters_vs_itr.txt).
import sys
import os
import statistics

run, first, last, step = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
SOLV = {"WAT", "HOH"}
IONS = {"NH4", "ACY", "Na+", "Cl-", "K+"}


def flag(lines, name):
    i = [l.strip() for l in lines].index("%FLAG " + name) + 2
    out = []
    while not lines[i].startswith("%"):
        out.append(lines[i])
        i += 1
    return out


top = open(os.path.join(run, "padded.parm7")).read().splitlines()
resname = [l[i:i + 4].strip() for l in flag(top, "RESIDUE_LABEL") for i in range(0, len(l), 4)]
ptr = [int(l[i:i + 8]) for l in flag(top, "RESIDUE_POINTER") for i in range(0, len(l), 8)]
atname = [l[i:i + 4].strip() for l in flag(top, "ATOM_NAME") for i in range(0, len(l), 4)]
ptr.append(len(atname) + 1)
resatoms = [atname[ptr[r] - 1:ptr[r + 1] - 1] for r in range(len(resname))]


def cls(res, at):
    rn = resname[res - 1]
    if rn in SOLV:
        return "water O" if at == "O" else "water H/EP"
    if rn in IONS:
        return "NH4/acetate"
    if at[0] == "H" or (at[0].isdigit() and at[1:2] == "H"):
        return "protein H"
    if at == "CA":
        return "protein CA"
    if at in ("N", "C", "O", "OXT"):
        return "protein N/C/O"
    return "protein side chain"


order = ["protein CA", "protein N/C/O", "protein side chain", "protein H", "water O", "water H/EP", "NH4/acetate"]
tot = {c: 0 for c in order}
nprot = sum(1 for r in resname if r not in SOLV | IONS)
for r, rn in enumerate(resname):
    if rn in SOLV | IONS:
        continue
    for a in resatoms[r]:
        tot[cls(r + 1, a)] += 1
nw = {}
for l in open(os.path.join(run, "nwaters_vs_itr.txt")):
    f = l.split()
    if f and f[0].isdigit():
        nw[int(f[0])] = int(f[1])

per = {c: [] for c in order}   # per iteration: list of weights of the restrained atoms in that class
nit = 0
for itr in range(first, last + 1, step):
    fn = os.path.join(run, "amber_%d.in" % itr)
    if not os.path.exists(fn):
        continue
    nit += 1
    L = open(fn).read().splitlines()
    w = {c: [] for c in order}
    i = 0
    while i < len(L):
        if L[i].strip().startswith("Specific"):
            k = float(L[i + 1])
            at = L[i + 3].split()[0]
            j = i + 5
            while not L[j].startswith("END"):
                f = L[j].split()
                for res in range(int(f[1]), int(f[2]) + 1):
                    if at in resatoms[res - 1]:
                        w[cls(res, at)].append(k)
                j += 1
            i = j
        i += 1
    for c in order:
        per[c].append(w[c])

wn = [nw[i] for i in range(first, last + 1, step) if i in nw]
print("%s: %d iterations %d-%d; %d protein residues; waters in the system %s" % (
    os.path.basename(run.rstrip("/")), nit, first, last, nprot, "%.0f" % statistics.mean(wn) if wn else "?"))
print("  %-19s %9s %11s %9s %9s %9s %9s" % ("atoms", "in system", "restrained", "k>=0.1", "k>=1", "median k", "max k"))
for c in order:
    if not any(per[c]):
        continue
    allw = [x for it in per[c] for x in it]
    n = statistics.mean(len(it) for it in per[c])
    n01 = statistics.mean(sum(1 for x in it if x >= 0.1) for it in per[c])
    n1 = statistics.mean(sum(1 for x in it if x >= 1) for it in per[c])
    sysn = tot.get(c, 0) if c.startswith("protein") else "-"
    print("  %-19s %9s %11.0f %9.0f %9.0f %9.3g %9.3g" % (c, sysn, n, n01, n1, statistics.median(allw), max(allw)))
