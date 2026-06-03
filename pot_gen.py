from params import *
import numpy as np
from matplotlib import pyplot as plot
import matplotlib
import os


matplotlib.rcParams["figure.dpi"] = 300
os.environ["PATH"] += os.pathsep + """/usr/local/texlive/2023/bin/
                                      universal-darwin"""
matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rc("font", family="STIXGeneral")
matplotlib.rc("font", weight="ultralight")


def force(phi: function, r: float, atom_i=None, atom_j=None) -> float:
    dr = 10**(-10)

    if atom_i is not None:
        return -(phi(r+dr, atom_i, atom_j) - phi(r, atom_i, atom_j)) / dr
    else:
        return -(phi(r+dr) - phi(r)) / dr


def eps_ij(atom_i: str, atom_j: str) -> float:
    return eps_hp * ((lambdas[atom_i] * lambdas[atom_j])**alpha)**(1/2)


def phi_hp(r: float, atom_i: str, atom_j: str) -> float:
    if r < sigma:
        return eps_rep * (sigma / r)**8 - eps_ij(atom_i, atom_j) * \
               (4 / 3 * (sigma / r)**6 - 1 / 3)

    return (eps_rep - eps_ij(atom_i, atom_j)) * (sigma / r)**8


def switch(r: float) -> float:
    return 1 - (10 * r**3 * sigma_hb**2 - 15 * r**4 * sigma_hb + 6 * r**5) / \
           sigma_hb**5


def phi_hb(r: float) -> float:
    return eps_hb * ((sigma_hb / (r + sigma_hb))**12 - 2 * (sigma_hb / (r + sigma_hb))**6) * switch(r)


def write_zero():
    table = open("ZERO.table", "w")
    table.write("ZERO\n")
    table.write("N 1001\n")
    table.write("\n")

    n = 1
    table.write(str(n) + " " + str(10**(-9)) + " " + str(0) + " " + \
                str(0) + "\n")

    rs = np.arange(0.01, 10.00 + 0.01, 0.01)

    for r in rs:
        n += 1
        table.write(str(n) + " " + str(r) + " " + str(0) + " " + str(0) + "\n")

    table.close()


def write_phi_hp(atom_i: str, atom_j: str):
    table = open(atom_i + "-" + atom_j + ".table", "w")
    table.write("BB\n")
    table.write("N 1001\n")
    table.write("\n")

    n = 1

    table.write(str(n) + " " + str(10**(-9)) + " " + \
                str(phi_hp(10**(-9), atom_i, atom_j)) + " " + \
                str(force(phi_hp, 10**(-9), atom_i, atom_j)) + "\n")

    rs = np.arange(0.01, 10.00 + 0.01, 0.01)

    for r in rs:
        n += 1
        table.write(str(n) + " " + str(r) + " " + \
                    str(phi_hp(r, atom_i, atom_j)) + " " + \
                    str(force(phi_hp, r, atom_i, atom_j)) + "\n")

    table.close()


def write_phi_hb():
    table = open("OH.table", "w")
    table.write("HB\n")
    table.write("N 1001\n")
    table.write("\n")
    n = 1
    table.write(str(n) + " " + str(10**(-9)) + " " + str(phi_hb(10**(-9))) + \
                " " + str(force(phi_hb, 10**(-9))) + "\n")
    rs = np.arange(0.01, 10.00 + 0.01, 0.01)
    for r in rs:
        n += 1
        table.write(str(n) + " " + str(r) + " " + str(phi_hb(r)) + " " + \
                    str(force(phi_hb, r)) + "\n")
    table.close()


write_zero()
write_phi_hp("BB", "BB")
write_phi_hp("BB", "SC")
write_phi_hp("SC", "SC")
write_phi_hb()


def draw_pot(table_name):
    f = open(table_name + ".table", "r")
    lines = f.readlines()[100:500]
    f.close()
    xs, ys = [], []
    for line in lines:
        nums = line.split()
        xs.append(float(nums[1]))
        ys.append(float(nums[2]))
    plot.plot(xs, ys, color="black")
    plot.xticks([])
    plot.yticks([])
    plot.xlabel("r, Å", fontsize=16)
    plot.ylabel("U(r), kcal / mol", fontsize=16)
    plot.tight_layout()
    #plot.show()
    plot.savefig("fig/" + table_name + ".png")

draw_pot("BB-BB")
draw_pot("BB-SC")
draw_pot("SC-SC")
