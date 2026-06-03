from math import log as ln
from random import uniform


seq = "GVVVVV"
seq = "GGGGGG"

N = 0
b = 0.5
M = 20
cell_size = 100
atom_types = 5
len_CP = 6.0
len_PO = len_PH = 2.0
L = len(seq)
L_res = L - seq.count("G")
print("The peptide contains", L, "residues and", L_res, "side chains.")

charge = {
    1: 0, 2: -2/3, 3: +2/3, 4: 0, 5: 0
}

radius = {
    1: 2.0, 2: 0.5, 3: 0.5, 4: 0.5, 5: 1.0
}

bead_nums = {
    "P": 1,
    "O": 2,
    "H": 3,
    "C": 4,
    "R": 5
}

f = open("data", "w")

f.write(" Self-assembly\n")
f.write("\n")
N_atoms = 2 * L - 1 + L_res + 2 * (L - 1)
f.write(str(M * (2 * L - 1 + L_res + 2 * (L - 1))) + " atoms\n")
f.write(str(M * (2 * L - 2 + L_res + 2 * (L - 1))) + " bonds\n")
f.write(str(M * (L - 1 + L - 1 + 2 * L - 3 + L - 1)) + " angles\n")
f.write("0 dihedrals\n")
f.write("0 impropers\n\n")

f.write("5 atom types\n")
f.write("3 bond types\n")
f.write("4 angle types\n")
f.write("0 dihedral types\n")
f.write("0 improper types\n\n")

f.write("0 " + str(cell_size) + " xlo xhi\n")
f.write("0 " + str(cell_size) + " ylo yhi\n")
f.write("0 " + str(cell_size) + " zlo zhi\n\n")

f.write(" Masses\n\n")
f.write("1 100.0\n")
f.write("2 10.0\n")
f.write("3 10.0\n")
f.write("4 100.0\n")
f.write("5 100.0\n\n")

f.write(" Atoms\n\n")

# solvent
#for n in range(1, N+1):
#    x, y, z = \
#        uniform(0, cell_size), uniform(0, cell_size), uniform(0, cell_size)
#    f.write(str(n) + " 0 0 " + str(x) + " " + str(y) + " " + str(z) + "\n")

n = 0

bonds_dict, angles_dict = dict(), dict()
n_bond, n_angle = 0, 0


def add_bead(x_: float, y_: float, z_: float, m: int, bead_type: str,
             side_chain=True):
    global n, n_bond, n_angle

    if bead_type == "P":
        x_P, y_P, z_P = x_, y_, z_
        x_O, y_O, z_O = x_, y_, z_ + len_PO
        x_H, y_H, z_H = x_, y_, z_ - len_PH

        n += 1
        f.write(str(n) + " " + str(m) + " 1 " + str(charge[1]) + " " + \
                str(x_P) + " " + str(y_P) + " " + str(z_P) + "\n")

        n_P = n

        n += 1
        n_bond += 1
        bonds_dict[n_bond] = (n - 1, n, 1)
        f.write(str(n) + " " + str(m) + " 2 " + str(charge[2]) + " " + \
                str(x_O) + " " + str(y_O) + " " + str(z_O) + "\n")

        n_O = n

        n += 1
        n_bond += 1
        bonds_dict[n_bond] = (n - 2, n, 1)
        f.write(str(n) + " " + str(m) + " 3 " + str(charge[3]) + " " + \
                str(x_H) + " " + str(y_H) + " " + str(z_H) + "\n")

        n_angle += 1
        angles_dict[n_angle] = (n - 1, n - 2, n, 1)

        return n_P, n_O
    elif bead_type == "Ca":
        x_Ca, y_Ca, z_Ca = x_, y_, z_
        x_R, y_R, z_R = x_, y_, z_ + radius[bead_nums["C"]] + radius[bead_nums["R"]]

        n += 1
        f.write(str(n) + " " + str(m) + " 4 " + str(charge[4]) + " " + \
                str(x_Ca) + " " + str(y_Ca) + " " + str(z_Ca) + "\n")

        n_Ca = n

        if side_chain:
            n += 1
            n_bond += 1
            bonds_dict[n_bond] = (n - 1, n, 3)
            f.write(str(n) + " " + str(m) + " 5 " + str(charge[5]) + " " + \
                    str(x_R) + " " + str(y_R) + " " + str(z_R) + "\n")

        return n_Ca
    else:
        raise ValueError("Unknown polypeptide bead type.")


# Beads.

C_alphas, O_bbs, Os, Hs, H_bbs = [], [], [], [], []

n_C_to_backbone_beads = dict()

for m in range(1, M+1):
    x_, y_, z_ = \
        uniform(0, cell_size), uniform(0, cell_size), uniform(0, cell_size)

    # Place Ca and peptide bond (P) beads evenly Ca P Ca P Ca ...

    backbone_beads = []
    for l in range(L):
        if seq[l] != "G":
            side_chain = True
        else:
            side_chain = False
        n_C = add_bead(
            x_ + len_CP * 2 * l,
            y_, z_, m, "Ca", side_chain)
        backbone_beads.append(n_C)
        n_C_to_backbone_beads[n_C] = len(backbone_beads) - 1
        if l < L - 1:
            n_P, n_O = add_bead(
                x_ + (len_CP) * (2 * l + 1),
                y_, z_, m, "P")
            backbone_beads.append(n_P)

            n_angle += 1
            angles_dict[n_angle] = (n_C, n_P, n_O, 3)

    # Establish bonds between consequent backbone beads.

    for l in range(len(backbone_beads) - 1):
        n_bond += 1
        bonds_dict[n_bond] = (backbone_beads[l], backbone_beads[l+1], 2)
        if l > 0:
            n_angle += 1
            if l % 2:
                angles_dict[n_angle] = (backbone_beads[l-1],
                                        backbone_beads[l],
                                        backbone_beads[l+1],
                                        2)
            else:
                angles_dict[n_angle] = (backbone_beads[l-1],
                                        backbone_beads[l],
                                        backbone_beads[l+1],
                                        4)

# Additional angles for orthogonality of h-bonding pseudoatoms to C-alpha
# planes.

n_angle_add = n_angle
angles_add_dict = dict()
for n in range(1, n_angle + 1):
    if angles_dict[n][-1] == 3:
        if n_C_to_backbone_beads[angles_dict[n][0]] + 2 < len(backbone_beads):
            n_angle_add += 1
            angles_add_dict[n_angle_add] = (n_C_to_backbone_beads[angles_dict[n][0]] + 2, n_P, n_O, 3)

angles_dict |= angles_add_dict

# Bonds.

f.write("\n Bonds\n\n")

for n_bond in bonds_dict:
    f.write(str(n_bond) + " " + str(bonds_dict[n_bond][2]) + " " + \
                                str(bonds_dict[n_bond][0]) + " " + \
                                str(bonds_dict[n_bond][1]) + "\n")

# Angles.

f.write("\n Angles\n\n")

for n_angle in angles_dict:
    angle = angles_dict[n_angle]
    f.write(str(n_angle) + " " + str(angle[3]) + " " + str(angle[0]) + " " + \
                                                       str(angle[1]) + " " + \
                                                       str(angle[2]) + "\n")

f.close()

# Excluded volumes.

f = open("settings_auto", "w")

f.write("pair_coeff * * morse 0.0 1.0 4.0 3.3068528194400546\n")
for i in range(1, atom_types+1):
    for j in range(i, atom_types+1):
        #f.write("pair_coeff " + str(i) + " " + str(j) + \
        #        " lj/cut 1.0 " + str(radius[i] + radius[j]) + " " + \
        #        str(radius[i] + radius[j]) + "\n")
                #str(round((radius[i] + radius[j]) / 1.12, 2)) + "\n")
        if i == j == 1 or i == j == 5:
            f.write("pair_coeff " + str(i) + " " + str(j) + \
                    " morse 10.0 1.0 " + str(radius[i] + radius[j] + ln(2) / 1.0) + " " + \
                    str(radius[i] + radius[j]) + "\n")

f.close()
