#!/usr/bin/env python3
import math
import random
import os

# ---------------------------
# Global parameters (if needed)
# ---------------------------
seed = None
mcsteps = None
d_max_mc = None

# ---------------------------
# Utility Functions
# ---------------------------
def fRand(fMin, fMax):
    """Return a random float between fMin and fMax."""
    return random.uniform(fMin, fMax)

# ---------------------------
# File output functions
# ---------------------------
def write_simCard(config, gamma, omegacc, omegacw, alpha, nsteps, ninfo, Lx, Ly, Lz,
                    nsubsteps, bc, margin, relax_time, nphases, mu, lambda_, kappa,
                    rad, xi, wallThich, wallKappa, SPOL, DPOL, JPOL, KPOL,
                    zetaS, zetaQ, SNEM, KNEM, JNEM, WNEM, count):
    """
    Write the simulation card file with the given parameters.
    """
    filename = f"simCard_{count}.dat"
    with open(filename, "w") as f:
        f.write("# Sample runcard\n")
        f.write(f"config = {config}\n")
        f.write(f"nsteps = {nsteps}\n")
        f.write(f"ninfo = {ninfo}\n")
        f.write(f"LX = {Lx}\n")
        f.write(f"LY = {Ly}\n")
        f.write(f"LZ = {Lz}\n")
        f.write(f"nsubsteps = {nsubsteps}\n")
        f.write(f"bc = {bc}\n")
        f.write(f"margin = {margin}\n")
        f.write(f"relax-time = {relax_time}\n")
        f.write(f"nphases = {nphases}\n")
        f.write(f"gamma = {gamma}\n")
        f.write(f"mu = {mu}\n")
        f.write(f"lambda = {lambda_}\n")
        f.write(f"kappa_cc = {kappa}\n")
        f.write(f"R = {rad}\n")
        f.write(f"xi = {xi}\n")
        f.write(f"omega_cc = {omegacc}\n")
        f.write(f"wall-thickness = {wallThich}\n")
        f.write(f"kappa_cs = {wallKappa}\n")
        f.write(f"omega_cs = {omegacw}\n")
        f.write(f"alpha = {alpha}\n")
        f.write(f"S-pol = {SPOL}\n")
        f.write(f"D-pol = {DPOL}\n")
        f.write(f"J-pol = {JPOL}\n")
        f.write(f"K-pol = {KPOL}\n")
        f.write(f"zetaS = {zetaS}\n")
        f.write(f"zetaQ = {zetaQ}\n")
        f.write(f"S-nem = {SNEM}\n")
        f.write(f"K-nem = {KNEM}\n")
        f.write(f"J-nem = {JNEM}\n")
        f.write(f"W-nem = {WNEM}\n")

def write_lattice(name, npart, xcf, ycf, zcf, zcoor):
    """
    Write the lattice positions to a file.
    Each line contains: x y zcoor (z coordinate is constant).
    """
    with open(name, "w") as f:
        for j in range(npart):
            f.write(f"{xcf[j]} {ycf[j]} {zcoor}\n")

def write_posfile_mix_perc(np, xcf, ycf, zcf, zcoor, count):
    """
    Write position file for mixed percentages.
    Always writes zcoor as the z-coordinate.
    """
    filename = "input_str.dat"  # Overwrites each time, similar to "w+"
    with open(filename, "w") as f:
        for j in range(np):
            f.write(f"{xcf[j]} {ycf[j]} {zcoor}\n")

def write_summary(count, gam1, gam2, zetas1, zetas2, omegacc1, omegacc2,
                  omegacw1, omegacw2, alpha, xi, zetaQ):
    """
    Append simulation parameter summary to file.
    """
    filename = "simulation_parameter_summary.dat"
    with open(filename, "a") as f:
        f.write(f"{count} {gam1} {gam2} {zetas1} {zetas2} {omegacc1} {omegacc2} "
                f"{omegacw1} {omegacw2} {alpha} {xi} {zetaQ}\n")

def Export(filename, tmp):
    """
    Export list of numbers to file, one per line.
    """
    with open(filename, "w") as f:
        for val in tmp:
            f.write(f"{val}\n")

# ---------------------------
# Lattice Initialization Functions
# ---------------------------
def init_square_lattice(l, npart, nx, ny, nz):
    """
    Create a square lattice of cell centers.
    Each cell center is computed as round(l*x + l/2) for x, and similarly for y and z.
    Returns lists xcf, ycf, zcf.
    """
    xcf = []
    ycf = []
    zcf = []
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                xc = round(l * x + l/2)
                yc = round(l * y + l/2)
                zc = round(l * z + l/2)
                xcf.append(xc)
                ycf.append(yc)
                zcf.append(zc)
    return xcf, ycf, zcf

def init_triangular_lattice(nx, ny, domainWidth, domainHeight, zcoor, R0):
    """
    Create a triangular lattice of cell centers.
    For even rows: x = i*dx, for odd rows: x = i*dx + 0.75*dx.
    Returns lists xcf, ycf, zcf.
    """
    dx = domainWidth / nx
    dy = domainHeight / ny
    dx2 = dx * 0.75  # offset for alternate rows
    xcf = []
    ycf = []
    zcf = []
    num = 0
    for j in range(ny):
        for i in range(nx):
            x = i * dx
            if j % 2 == 1:
                x += dx2
            y = j * dy
            xcf.append(x)
            ycf.append(y + R0)
            zcf.append(zcoor)
            num += 1
    print("num:", num)
    return xcf, ycf, zcf

def disorder_initial(l, nx, ny, nz):
    """
    Create a disordered lattice: for each grid point compute its center and return lists.
    """
    xcf = []
    ycf = []
    zcf = []
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                xc = round(l * x + l/2)
                yc = round(l * y + l/2)
                zc = round(l * z + l/2)
                xcf.append(xc)
                ycf.append(yc)
                zcf.append(zc)
    return xcf, ycf, zcf

def disorder_mc(d_max, xcf, ycf, zcf):
    """
    Apply a Monte Carlo displacement to lattice positions.
    """
    npart = len(xcf)
    for i in range(npart):
        exp1 = random.randint(0, 1)
        exp2 = random.randint(0, 1)
        exp3 = random.randint(0, 1)
        p1 = fRand(0, d_max)
        p2 = fRand(0, d_max)
        p3 = fRand(0, d_max)
        xcf[i] = xcf[i] + ((-1) ** exp1) * p1
        ycf[i] = ycf[i] + ((-1) ** exp2) * p2
        zcf[i] = zcf[i] + ((-1) ** exp3) * p3

# ---------------------------
# Other Utility Functions
# ---------------------------
def find_min(vect):
    """Return the smallest element in the list."""
    return min(vect)

def compute_all_dist(index, ix, iy, iz, txc, tyc, tzc):
    """
    Compute the minimum distance from point (ix,iy,iz) to all other points in txc,tyc,tzc.
    """
    distances = []
    npart = len(txc)
    for j in range(npart):
        if j != index:
            d = math.sqrt((txc[j] - ix)**2 + (tyc[j] - iy)**2 + (tzc[j] - iz)**2)
            distances.append(d)
    if distances:
        return find_min(distances)
    else:
        return None

# ---------------------------
# Main Routine
# ---------------------------
def main():
    # Simulation parameters
    nsteps = 1000
    ninfo = 10
    Lz = int(40.0)

    ncells = 9   # total number of cells (e.g., 9 => 3x3 grid)
    R0 = 8.0
    rad = R0
    lbox = 2.0 * R0
    nsqrt = int(math.sqrt(ncells))
    Lx = int(nsqrt * lbox)
    Ly = int(nsqrt * lbox)

    nsubsteps = 10
    bc = 2
    margin = int(18.0)
    relax_time = 100
    wall_thickness = 7.0
    zcoor = wall_thickness + R0/2.0

    gamma = 0.008
    omegacc = 0.0008
    omegacw = 0.002
    alpha = 0.05
    xi = 1.0
    zetaS = 0.0
    zetaQ = 0.0

    mu = 45.0
    lambda_ = 3.0
    kappa = 0.5

    wall_kappa = 0.2
    SPOL = 1.0
    SNEM = 0.0
    KNEM = 0.0
    JNEM = 0.0
    WNEM = 0.0
    JPOL = 0.001
    KPOL = 0.001
    DPOL = 0.001

    global d_max_mc
    d_max_mc = 1.0

    npx = nsqrt
    npy = nsqrt
    npz = 1
    nphases = int(ncells)

    print(f"Number of cells are: {nphases} with rad: {rad}")

    # Allocate lattice arrays (using lists)
    total_points = npx * npy * npz
    # Initialize with zeros (or you can later overwrite these values)
    xcf = [0.0] * total_points
    ycf = [0.0] * total_points
    zcf = [0.0] * total_points

    # Initialize a square lattice
    xcf, ycf, zcf = init_square_lattice(lbox, nphases, npx, npy, npz)
    # Optionally, you can apply disorder:
    # disorder_mc(d_max_mc, xcf, ycf, zcf)
    # For debugging, you might write out the lattice to a file
    # write_lattice("lattice.dat", nphases, xcf, ycf, zcf, zcoor)

    config = "input const"

    # Define mixing parameter arrays (each with one element)
    zetaS_mix_a = [0.0]  # wt 
    zetaS_mix_b = [0.0]  # ecad 
    gamma_mix_a = [0.008]
    gamma_mix_b = [0.008]
    omegacc_mix_a = [0.0008]  # wt 
    omegacc_mix_b = [0.0008]  # ecad 
    omegacw_mix_a = [0.0020]  # wt 
    omegacw_mix_b = [0.0020]  # ecad 
    zetaQ_mix_a = [0.0]  # wt 
    zetaQ_mix_b = [0.0]  # ecad 
    alpha_mix = [alpha]
    xi_mix = [xi]
    kappa_mix = [kappa]
    mu_mix = [mu]
    rad_mix = [rad]

    count = 1

    # Nested loops (since each array has one element, this loop runs once)
    for i in range(len(gamma_mix_a)):
        for j in range(len(omegacc_mix_a)):
            for k in range(len(omegacw_mix_a)):
                for m in range(len(alpha_mix)):
                    for ii in range(len(zetaS_mix_a)):
                        for jj in range(len(zetaQ_mix_a)):
                            for kk in range(len(xi_mix)):
                                for mm in range(len(kappa_mix)):
                                    for ll in range(len(mu_mix)):
                                        for hh in range(len(rad_mix)):
                                            GAMMA_A = gamma_mix_a[i]
                                            GAMMA_B = gamma_mix_b[i]

                                            ZETAS_A = zetaS_mix_a[ii]
                                            ZETAS_B = zetaS_mix_b[ii]

                                            OMEGACC_A = omegacc_mix_a[j]
                                            OMEGACC_B = omegacc_mix_b[j]

                                            OMEGACW_A = omegacw_mix_a[k]
                                            OMEGACW_B = omegacw_mix_b[k]

                                            ZETAQ_A = zetaQ_mix_a[jj]
                                            ZETAQ_B = zetaQ_mix_b[jj]

                                            ALPHA = alpha_mix[m]
                                            XI = xi_mix[kk]
                                            KAPPA = kappa_mix[mm]
                                            MU = mu_mix[ll]
                                            RAD = rad_mix[hh]

                                            # Write position file (always writes "input_str.dat")
                                            write_posfile_mix_perc(int(ncells), xcf, ycf, zcf, zcoor, count)

                                            # Write simulation card file
                                            write_simCard(config, gamma, omegacc, omegacw, alpha,
                                                          nsteps, ninfo, Lx, Ly, Lz, nsubsteps, bc,
                                                          margin, relax_time, int(ncells), mu, lambda_,
                                                          kappa, rad, xi, wall_thickness, wall_kappa,
                                                          SPOL, DPOL, JPOL, KPOL, zetaS, zetaQ,
                                                          SNEM, KNEM, JNEM, WNEM, count)

                                            # You can compute ratios if needed:
                                            ratio_a = OMEGACC_A / OMEGACC_B if OMEGACC_B != 0 else None
                                            ratio_b = OMEGACW_A / OMEGACW_B if OMEGACW_B != 0 else None

                                            # Optionally, write summary (commented out in original)
                                            # write_summary(count, GAMMA_A, GAMMA_B, ZETAS_A, ZETAS_B,
                                            #               OMEGACC_A, OMEGACC_B, OMEGACW_A, OMEGACW_B,
                                            #               ALPHA, XI, ZETAQ_A)
                                            count += 1

    print("Computation done.")

if __name__ == "__main__":
    main()

