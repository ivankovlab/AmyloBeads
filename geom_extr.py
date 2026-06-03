import numpy as np
from Bio.PDB import PDBParser


amyloid = PDBParser().get_structure("VILLIS", "VILLIS.pdb")
chains = ["B", "A", "C", "E", "G"]

Cs, Os, Ns, Hs = dict(), dict(), dict(), dict()
for model in amyloid:
    for chain in model:
        if chain._id in chains:
            Cs[chain._id], Os[chain._id], Ns[chain._id], Hs[chain._id] = list(), list(), list(), list()
            for atom in chain.get_atoms():
                if atom.get_name() == "O":
                    Os[chain._id].append(atom.get_coord())
                elif atom.get_name() == "N":
                    Ns[chain._id].append(atom.get_coord())
                elif atom.get_name()[0] == "H":
                    Hs[chain._id].append(atom.get_coord())
                elif atom.get_name() == "CA":
                    Cs[chain._id].append([float(atom.get_coord()[0]),
                                          float(atom.get_coord()[1]),
                                          float(atom.get_coord()[2])])


def get_dist(a, b):
    return((a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2)**(1/2)


def get_angle(a, b, c):
    a, b, c = np.array(a, dtype=float), np.array(b, dtype=float), np.array(c, dtype=float)
    u, v = a - b, c - b
    return np.acos(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v)) * 180 / np.pi


for _, chain in Cs.items():
    print(chain, ",")


min_dist = None
for C_1 in Cs["A"]:
    for C_2 in Cs["B"]:
        dist = get_dist(C_1, C_2)
        if min_dist is None or dist < min_dist:
            min_dist = dist

print("C-C between strands", round(min_dist, 2))

min_dist = None
for C_1 in Cs["A"]:
    for C_2 in Cs["C"]:
        dist = get_dist(C_1, C_2)
        if min_dist is None or dist < min_dist:
            min_dist = dist

print("C-C", round(min_dist, 2))

min_dist = None
for O in Os["A"]:
    for N in Ns["C"]:
        dist = get_dist(O, N)
        if min_dist is None or dist < min_dist:
            min_dist = dist

print("N-O", round(min_dist, 2))

min_dist = None
for O in Os["A"]:
    for H in Hs["C"]:
        dist = get_dist(O, H)
        if min_dist is None or dist < min_dist:
            min_dist = dist

print("H-O", round(min_dist, 2))
