from params import *
from math import log as ln
from random import uniform
import numpy as np


# Parameters of the model.
seq = "G" * 20
cell_size = 200#100
charge = {1: 0, 2: 0, 3: -2/3, 4: 0, 5: +2/3, 6: 0}

L = len(seq)

chains = [
[[5.866000175476074, 6.76800012588501, 4.058000087738037], [8.673999786376953, 7.085999965667725, 6.586999893188477], [12.005999565124512, 6.854000091552734, 4.75], [15.633999824523926, 7.206999778747559, 5.796000003814697], [18.54800033569336, 7.179999828338623, 3.321000099182129], [22.20800018310547, 7.310999870300293, 4.311999797821045]] ,
[[5.866000175476074, 1.9299999475479126, 4.058000087738037], [8.673999786376953, 2.247999906539917, 6.586999893188477], [12.005999565124512, 2.0160000324249268, 4.75], [15.633999824523926, 2.36899995803833, 5.796000003814697], [18.54800033569336, 2.3420000076293945, 3.321000099182129], [22.20800018310547, 2.4730000495910645, 4.311999797821045]] ,
[[5.866000175476074, -2.9079999923706055, 4.058000087738037], [8.673999786376953, -2.5899999141693115, 6.586999893188477], [12.005999565124512, -2.822000026702881, 4.75], [15.633999824523926, -2.4690001010894775, 5.796000003814697], [18.54800033569336, -2.496000051498413, 3.321000099182129], [22.20800018310547, -2.365000009536743, 4.311999797821045]] ,
[[5.866000175476074, -7.745999813079834, 4.058000087738037], [8.673999786376953, -7.427999973297119, 6.586999893188477], [12.005999565124512, -7.659999847412109, 4.75], [15.633999824523926, -7.307000160217285, 5.796000003814697], [18.54800033569336, -7.334000110626221, 3.321000099182129], [22.20800018310547, -7.203000068664551, 4.311999797821045]]
]

chains = None

N_beta, N_anch = 12, 4
origin = (30, 30, 30)
M_free = 1

if chains is not None:
    M = len(chains)
    L = len(chains[0])

if (N_beta is not None) and (N_anch is not None):
    M = N_beta * N_anch + M_free
    seed = dict()
    m = 0
    for n_anch in range(N_anch):
        for n_beta in range(N_beta):
            m += 1
            seed[m] = (origin[0],
                       origin[1] + n_beta * d_beta,
                       origin[2] + n_anch * d_anch)

atom_types = 6

atoms_dict, bonds_dict, angles_dict, dihedrals_dict = dict(), dict(), dict(), dict()
n_atom, n_bond, n_angle, n_dihedral = 0, 0, 0, 0

# Add polypeptides.

for m in range(1, M+1):
    C_alphas, O_bbs, Os, H_bbs, Hs, Rs = [], [], [], [], [], []

    if (N_beta is not None) and (N_anch is not None):
        # If in a filament, define position by the lattice rule.
        if m <= N_beta * N_anch:
            x_, y_, z_ = seed[m]
        # Free monomers in vicinity of the tip.
        else:
            x_ = uniform(35, 40)
            y_ = uniform(20, 30)
            z_ = uniform(30, 40)
    else:
        x_, y_, z_ = \
            uniform(0, cell_size), uniform(0, cell_size), uniform(0, cell_size)

    dir = -1

    for l in range(L):
        if chains is None:
            x_dir, z_dir = l_bb * np.sin(60 / 180 * np.pi), -l_bb * 0.5 * dir
            N_ter_x = x_ + l_bb * np.sin(60 / 180 * np.pi) * l
            N_ter_y = y_
            if dir < 0:
                N_ter_z = z_
            else:
                N_ter_z = z_ + l_bb * 0.5
        else:
            N_ter_x = chains[m-1][l][0] + 50
            N_ter_y = chains[m-1][l][1] + 50
            N_ter_z = chains[m-1][l][2] + 50
            if l < L - 1:
                x_dir = (chains[m-1][l+1][0] - chains[m-1][l][0])
                z_dir = (chains[m-1][l+1][2] - chains[m-1][l][2])

        dir = -dir

        n_atom += 1
        atoms_dict[n_atom] = (m, 1, charge[1], N_ter_x, N_ter_y, N_ter_z)
        C_alphas.append(n_atom)

        # Add peptide bond, if the C-alpha atom is not the last one.
        if l < L - 1:
            n_atom += 1
            atoms_dict[n_atom] = (m, 2, charge[2],
                N_ter_x + x_dir / 3, N_ter_y, N_ter_z + z_dir / 3)
            O_bbs.append(n_atom)

            n_atom += 1
            atoms_dict[n_atom] = (m, 3, charge[3],
                N_ter_x + x_dir / 3, N_ter_y + dir * l_hb, N_ter_z + z_dir / 3)
            Os.append(n_atom)

            n_atom += 1
            atoms_dict[n_atom] = (m, 4, charge[4],
                N_ter_x + x_dir * 2 / 3, N_ter_y, N_ter_z + z_dir * 2 / 3)
            H_bbs.append(n_atom)

            n_atom += 1
            atoms_dict[n_atom] = (m, 5, charge[5],
                N_ter_x + x_dir * 2 / 3, N_ter_y - dir * l_hb, N_ter_z + z_dir * 2 / 3)
            Hs.append(n_atom)

            if 0 < l < L - 1:
                n_atom += 1
                atoms_dict[n_atom] = (m, 6, charge[6],
                    N_ter_x, N_ter_y, N_ter_z - dir * l_sc)
                Rs.append(n_atom)

    for i in range(len(C_alphas) - 1):
        n_bond += 1
        bonds_dict[n_bond] = (C_alphas[i], C_alphas[i+1], 1)

        n_bond += 1
        bonds_dict[n_bond] = (C_alphas[i], O_bbs[i], 2)

        n_bond += 1
        bonds_dict[n_bond] = (H_bbs[i], C_alphas[i+1], 2)

        n_bond += 1
        bonds_dict[n_bond] = (O_bbs[i], Os[i], 3)

        n_bond += 1
        bonds_dict[n_bond] = (H_bbs[i], Hs[i], 3)

        if 0 < i < len(C_alphas) - 1:
            n_bond += 1
            bonds_dict[n_bond] = (C_alphas[i], Rs[i-1], 4)

        n_angle += 1
        angles_dict[n_angle] = (C_alphas[i], O_bbs[i], C_alphas[i+1], 1)

        n_angle += 1
        angles_dict[n_angle] = (C_alphas[i], H_bbs[i], C_alphas[i+1], 1)

        n_angle += 1
        angles_dict[n_angle] = (C_alphas[i], O_bbs[i], Os[i], 2)

        n_angle += 1
        angles_dict[n_angle] = (C_alphas[i], H_bbs[i], Hs[i], 2)

        if 0 < i < len(C_alphas) - 1:
            n_angle += 1
            angles_dict[n_angle] = (C_alphas[i-1], C_alphas[i], Rs[i-1], 3)
            n_angle += 1
            angles_dict[n_angle] = (C_alphas[i+1], C_alphas[i], Rs[i-1], 3)

        if i < L - 2:
            n_angle += 1
            angles_dict[n_angle] = (Os[i], O_bbs[i], C_alphas[i+2], 2)

            n_angle += 1
            angles_dict[n_angle] = (Hs[i], H_bbs[i], C_alphas[i+2], 2)

            n_angle += 1
            angles_dict[n_angle] = (C_alphas[i], C_alphas[i+1], C_alphas[i+2], 3)

        n_dihedral += 1
        dihedrals_dict[n_dihedral] = (Os[i], O_bbs[i], H_bbs[i], Hs[i], 1)

f = open("data", "w")

f.write(" Self-assembly\n")
f.write("\n")
f.write(str(len(atoms_dict)) + " atoms\n")
f.write(str(len(bonds_dict)) + " bonds\n")
f.write(str(len(angles_dict)) + " angles\n")
f.write(str(len(dihedrals_dict)) + " dihedrals\n")
f.write("0 impropers\n\n")

f.write("6 atom types\n")
f.write("4 bond types\n")
f.write("3 angle types\n")
f.write("1 dihedral types\n")
f.write("0 improper types\n\n")

f.write("0 " + str(cell_size) + " xlo xhi\n")
f.write("0 " + str(cell_size) + " ylo yhi\n")
f.write("0 " + str(cell_size) + " zlo zhi\n\n")

f.write(" Masses\n\n")
f.write("1 100.0\n")
f.write("2 10.0\n")
f.write("3 10.0\n")
f.write("4 10.0\n")
f.write("5 10.0\n")
f.write("6 100.0\n\n")

# Write atoms.

f.write(" Atoms\n\n")

for n_atom in atoms_dict:
    atom = atoms_dict[n_atom]
    f.write(str(n_atom) + " " + str(atom[0]) + " " + str(atom[1]) + " " + \
                                str(atom[2]) + " " + str(atom[3]) + " " + \
                                str(atom[4]) + " " + str(atom[5]) + "\n")

# Write bonds.

f.write("\n Bonds\n\n")

for n_bond in bonds_dict:
    f.write(str(n_bond) + " " + str(bonds_dict[n_bond][2]) + " " + \
                                str(bonds_dict[n_bond][0]) + " " + \
                                str(bonds_dict[n_bond][1]) + "\n")

# Write angles.

f.write("\n Angles\n\n")

for n_angle in angles_dict:
    angle = angles_dict[n_angle]
    f.write(str(n_angle) + " " + str(angle[3]) + " " + str(angle[0]) + " " + \
                                                       str(angle[1]) + " " + \
                                                       str(angle[2]) + "\n")

# Write dihedrals.

f.write("\n Dihedrals\n\n")

for n_dihedral in dihedrals_dict:
    dihedral = dihedrals_dict[n_dihedral]
    f.write(str(n_dihedral) + " " + str(dihedral[4]) + " " +
                                    str(dihedral[0]) + " " + \
                                    str(dihedral[1]) + " " + \
                                    str(dihedral[2]) + " " + \
                                    str(dihedral[3]) + "\n")

f.close()

# Set bonds parameters depending on backbone and side chain geometry.

f = open("settings_auto", "w")
f.write("# distance between consequent alpha-carbons\n")
f.write("bond_coeff 1 100 " + str(l_bb) + "\n")
f.write("# one-third\n")
f.write("bond_coeff 2 100 " + str(l_bb / 3) + "\n")
f.write("# radius of the pseudoatom\n")
f.write("bond_coeff 3 100 " + str(l_hb) + "\n")
f.write("# distance from Ca to side chain\n")
f.write("bond_coeff 4 100 " + str(l_sc) + "\n")
f.close()
