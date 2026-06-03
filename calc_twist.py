import numpy as np
import numpy


def get_L_vec(chains, m: int) -> numpy.array:
    """Calculate vector connecting N- and C-terminal beads of a polypeptide
       chain.

    Parameters
    ----------
        m : int
            Number of the chain.

    Returns
    -------
        Numpy array of three floats.
    """

    N_ter, C_ter = chains[m][0], chains[m][-1]
    N_ter_x, N_ter_y, N_ter_z = N_ter[1], N_ter[2], N_ter[3]
    C_ter_x, C_ter_y, C_ter_z = C_ter[1], C_ter[2], C_ter[3]
    #print(N_ter_x, N_ter_y, N_ter_z)
    #print(C_ter_x, C_ter_y, C_ter_z)
    return np.array([C_ter_x - N_ter_x, C_ter_y - N_ter_y, C_ter_z - N_ter_z],
                    dtype=float)


def get_bb_vec(chains, m: int, a: int) -> numpy.array:
    bb_bead, bb_bead_ = chains[m][a], chains[m][a+5]
    x, y, z = bb_bead[1], bb_bead[2], bb_bead[3]
    x_, y_, z_ = bb_bead_[1], bb_bead_[2], bb_bead_[3]
    return np.array([x_ - x, y_ - y, z_ - z], dtype=float)


def get_local_twist(chains, m: int, m_: int) -> float:
    """Calculate twist between two polypeptide chains in a filament.

    Parameters
    ----------
        m : int
            Number of the first chain.
        m_ : int
            Number of the second chain.

    Returns
    -------
        One float number.
    """

    L = 0
    for atom in chains[m]:
        if atom[4] == 1:
            L += 1

    #twist_sum = 0

    #for a in range(0, L-1):
    #    v, v_ = get_bb_vec(chains, m, a * 5), get_bb_vec(chains, m_, a * 5)
    #    cos_val = np.dot(v, v_) / np.linalg.norm(v) / np.linalg.norm(v_)
    #    cos_val = max(-1, cos_val)
    #    cos_val = min(1, cos_val)
    #    twist_sum += np.acos(cos_val) * 180 / np.pi

    #return twist_sum / (L - 1)

    v, v_ = get_L_vec(chains, m), get_L_vec(chains, m_)
    cos_val = np.dot(v, v_) / np.linalg.norm(v) / np.linalg.norm(v_)
    cos_val = max(-1, cos_val)
    cos_val = min(1, cos_val)

    return np.arccos(cos_val) * 180 / np.pi


def get_twist(data_file: str) -> float:
    # Read all lines from LAMMPS data file.
    f = open(data_file, "r")
    lines = f.readlines()
    f.close()

    chains = dict()

    atoms_flag = False

    for line in lines:
        tokens = line.split()

        if len(tokens) == 0: continue

        # Read atom's coordinates.
        if atoms_flag and len(tokens) >= 7:
            a, m, at = int(tokens[0]), int(tokens[1]), int(tokens[2])
            x, y, z = float(tokens[4]), float(tokens[5]), float(tokens[6])
            if m in chains:
                chains[m].append((a, x, y, z, at))
            else:
                chains[m] = [(a, x, y, z, at)]

        # Start reading "Atoms" section.
        if tokens[0] == "Atoms":
            atoms_flag = True

        # Finish reading "Atoms" section.
        if tokens[0] == "Bonds":
            atoms_flag = False

    # Order atoms in each chain by atom number.
    for m in chains:
        chains[m] = sorted(chains[m])

    # Calculate the local twists and their average.
    local_twists = []

    # The first strand.
    for m in range(2, 10+1):
        local_twist = get_local_twist(chains, m, m+1)
        local_twists.append(local_twist)
        #print(m, local_twist)

    # The second strand.
    for m in range(14, 22+1):
        local_twist = get_local_twist(chains, m, m+1)
        local_twists.append(local_twist)
        #print(m, local_twist)

    local_twists = np.array(local_twists, dtype=float)
    #print(round(np.mean(local_twists), 2))

    return np.mean(local_twists)


#for t in [0, 250, 500, 750, 1000]:
#    data_file = "res/" + str(t) + ".data"
#    print(get_twist(data_file))

print("        dir LJ LJ",
      round(get_twist("res_traj/dir_LJ_LJ.data") / 10, 2))
print("     dir Morse LJ",
      round(get_twist("res_traj/dir_Morse_LJ.data") / 10, 2))
print("     dir LJ Morse",
      round(get_twist("res_traj/dir_LJ_Morse.data") / 10, 2))
print("  dir Morse Morse",
      round(get_twist("res_traj/dir_Morse_Morse.data") / 10, 2))
print("      undir LJ LJ",
      round(get_twist("res_traj/undir_LJ_LJ.data") / 10, 2))
print("   undir LJ Morse",
      round(get_twist("res_traj/undir_LJ_Morse.data") / 10, 2))
print("   undir Morse LJ",
      round(get_twist("res_traj/undir_Morse_LJ.data") / 10, 2))
print("undir Morse Morse",
      round(get_twist("res_traj/undir_Morse_Morse.data") / 10, 2))
