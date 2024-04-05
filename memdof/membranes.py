import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.vectors import calc_dihedral
from mtypes import *
from scipy.signal import find_peaks, savgol_filter
from scipy.stats import normaltest

def gaussian(x, mu, sig):
    return 1./(np.sqrt(2.*np.pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2)

def parse_PTB(
    path: str, id: str = "X", ignore_hydrogen: bool = False, quiet: bool = True
) -> tuple[Structure, tuple[float, float, float]]:
    P = PDBParser(QUIET=quiet)

    structure: Structure = P.get_structure(id, path)
    pbc: tuple[float, float, float] = (0, 0, 0)

    # Read file to check for periodic boundary conditions
    with open(path, "r") as f:
        for line in f.readlines():
            if line.startswith("CRYST1"):
                pbc = [float(x) for x in line.split()[1:4]]

    return structure, pbc


if __name__ == "__main__":
    structure, pbc = parse_PTB("tests/assets/1.pdb", "DOPC", quiet=True)

    STATISTIC_THRESHOLDS = [115, 115, 105] # Bonds, angles, dihedrals
    MAX_RELATIVE_STDS = [0.01, 0.05, 0.05] # Bonds, angles, dihedrals

    import matplotlib.pyplot as plt
    print(plt.style.available)
    plt.style.use("seaborn-pastel")
    from topologies import parse_topology

    bond_info = parse_topology("tests/assets/DOPC.itp", ignore_hydrogen=True)

    norms = [[] for _ in bond_info.bonds]
    angles = [[] for _ in bond_info.angles]
    dihedrals = [[] for _ in bond_info.dihedrals]

    for chain in structure.get_chains():
        for residue in chain.get_residues():
            atoms = list(residue.get_atoms())

            # Iterate over all bonds
            for coordinate_index, coordinate in enumerate(bond_info.bonds):
                # Get atom names of the bond
                i_name, j_name = (
                    bond_info.atoms[coordinate["i"]]["atom"],
                    bond_info.atoms[coordinate["j"]]["atom"],
                )

                # Get atom objects
                i_atom, j_atom = (
                    [a for a in atoms if a.get_name() == i_name][0],
                    [a for a in atoms if a.get_name() == j_name][0],
                )

                # Calculate the vector between the atoms
                vec = i_atom.get_vector() - j_atom.get_vector()
                vec = np.array([vec[0], vec[1], vec[2]])

                # Apply periodic boundary conditions
                if np.linalg.norm(vec) > np.min(pbc):
                    vec = vec - pbc * np.round(vec / pbc)

                # Append the norm of the vector to the list
                norms[coordinate_index].append(np.linalg.norm(vec))

            # Iterate over all angles
            for angle_index, angle in enumerate(bond_info.angles):
                # Get atom names of the bond
                i_name, j_name, k_name = (
                    bond_info.atoms[angle["i"]]["atom"],
                    bond_info.atoms[angle["j"]]["atom"],
                    bond_info.atoms[angle["k"]]["atom"],
                )

                # Get atom objects
                i_atom, j_atom, k_atom = (
                    [a for a in atoms if a.get_name() == i_name][0],
                    [a for a in atoms if a.get_name() == j_name][0],
                    [a for a in atoms if a.get_name() == k_name][0],
                )

                # Calculate the two vectors between the atoms
                vec1 = i_atom.get_vector() - j_atom.get_vector()
                vec2 = k_atom.get_vector() - j_atom.get_vector()

                # Make sure the vectors are numpy arrays
                vec1 = np.array([vec1[0], vec1[1], vec1[2]])
                vec2 = np.array([vec2[0], vec2[1], vec2[2]])

                # Apply periodic boundary conditions
                if np.linalg.norm(vec1) > np.min(pbc):
                    vec1 = vec1 - pbc * np.round(vec1 / pbc)
                if np.linalg.norm(vec2) > np.min(pbc):
                    vec2 = vec2 - pbc * np.round(vec2 / pbc)

                # Calculate the angle between the vectors
                angle = np.arccos(
                    np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
                )

                # To degrees
                angle = np.degrees(angle)

                # Append the angle to the list
                angles[angle_index].append(angle)

            # Iterate over all dihedrals
            for dihedral_index, dihedral in enumerate(bond_info.dihedrals):
                # Get atom names of the bond
                i_name, j_name, k_name, l_name = (
                    bond_info.atoms[dihedral["i"]]["atom"],
                    bond_info.atoms[dihedral["j"]]["atom"],
                    bond_info.atoms[dihedral["k"]]["atom"],
                    bond_info.atoms[dihedral["l"]]["atom"],
                )

                # Get atom objects
                i_atom, j_atom, k_atom, l_atom = (
                    [a for a in atoms if a.get_name() == i_name][0],
                    [a for a in atoms if a.get_name() == j_name][0],
                    [a for a in atoms if a.get_name() == k_name][0],
                    [a for a in atoms if a.get_name() == l_name][0],
                )

                # Apply periodic boundary conditions
                vec1 = i_atom.get_vector() - j_atom.get_vector()
                vec2 = k_atom.get_vector() - j_atom.get_vector()
                vec3 = l_atom.get_vector() - k_atom.get_vector()

                # Make sure the vectors are numpy arrays
                vec1 = np.array([vec1[0], vec1[1], vec1[2]])
                vec2 = np.array([vec2[0], vec2[1], vec2[2]])
                vec3 = np.array([vec3[0], vec3[1], vec3[2]])

                if np.linalg.norm(vec1) > np.min(pbc):
                    vec1 = vec1 - pbc * np.round(vec1 / pbc)
                if np.linalg.norm(vec2) > np.min(pbc):
                    vec2 = vec2 - pbc * np.round(vec2 / pbc)
                if np.linalg.norm(vec3) > np.min(pbc):
                    vec3 = vec3 - pbc * np.round(vec3 / pbc)

                # Calculate the dihedral
                dihedral = np.arccos(
                    np.dot(
                        np.cross(vec1, vec2) / np.linalg.norm(np.cross(vec1, vec2)),
                        np.cross(vec2, vec3) / np.linalg.norm(np.cross(vec2, vec3)),
                    )
                    / (
                        np.linalg.norm(np.cross(vec1, vec2))
                        * np.linalg.norm(np.cross(vec2, vec3))
                    )
                )

                # To degrees
                dihedral = np.degrees(dihedral)

                # Append the dihedral to the list
                dihedrals[dihedral_index].append(dihedral)

    # Internal coodiantes
    coordinates = [norms, angles, dihedrals]
    fixed = [[], [], []]    # Fixed norms, angles and dihedrals
    freedom = [[], [], []]  # Free norms, angles and dihedrals
    coodinates_with_fixed_flag = [[], [], []]

    for k, coordinate_list in enumerate([bond_info.bonds, bond_info.angles, bond_info.dihedrals]):
        print("Starting analysis of ", "bonds" if k == 0 else "angles" if k == 1 else "dihedrals")

        STATISTIC_THRESHOLD = STATISTIC_THRESHOLDS[k]
        MAX_RELATIVE_STD = MAX_RELATIVE_STDS[k]

        for coordinate_index, coordinate in enumerate(coordinate_list):
            # Calculate the mean and standard deviation of the norms
            mean, std = np.mean(coordinates[k][coordinate_index]), np.std(coordinates[k][coordinate_index])

            # Bin the norms in 100 bins
            hist, bin_edges = np.histogram(coordinates[k][coordinate_index], bins=100)

            # Find value where the histogram is at its peak
            peak0 = np.argmax(hist)

            # Remove all norms that are less than 5% of the peaks occurrence
            hist[hist < 0.05 * hist[peak0]] = 0

            # Apply a Savitzky-Golay filter to the histogram
            hist = savgol_filter(hist, 21, 3)

            # Remove negative values
            hist[hist < 0] = 0

            # Calculate the normality of the data
            res = normaltest(coordinates[k][coordinate_index])
            normality = res.statistic

            # Decide if the data is normal (normal = fixed)
            is_fixed = normality < STATISTIC_THRESHOLD and std / mean < MAX_RELATIVE_STD

            # Find peaks for histogram
            peaks, dict = find_peaks(
                hist,
                threshold=[0, 99999],
                distance=20,
                width=10,
                height=[np.max(hist) * 0.05, 99999],
            )

            # Choose different colors for different coordiantes
            color = 'black'
            if k == 0:
                color = "purple"
            elif k == 1:
                color = "green"
            elif k == 2:
                color = "orange"

            plt.grid()

            # Plot the histogram and the peaks bonds of the original data
            plt.hist(coordinates[k][coordinate_index], bins=100, alpha=0.75, color=color, label="Bond length occurances")

            # Plot the smoothed histogram (same x-axis as the original data)
            plt.plot(bin_edges[1:], hist, color=color, label="Smoothed histogram", alpha=0.75)

            # Plot the theoretical normal distribution (continuous)
            x = np.linspace(np.min(coordinates[k][coordinate_index]), np.max(coordinates[k][coordinate_index]), 100)
            y = gaussian(x, mean, std)

            # Scale the normal distribution to the histogram
            y = y / np.max(y) * np.max(hist)

            # Choose color for different types of internal coordinates
            color = "b" if is_fixed else "r"
            opacity = 1 if is_fixed else 0.5
            dimension = "Å" if k == 0 else "°"

            # Plot normal distrubtion
            plt.plot(x, y, color=color, linestyle="dotted", label="Normal distribution", alpha=opacity, linewidth=1.5)

            # Plot the mean and standard deviation
            plt.axvline(mean, color="black", linestyle="--", alpha=0.5, linewidth=0.75, label=f"Mean $\mu={mean:.3f}{dimension}$")
            plt.axvline(mean + std, color="black", linestyle="--", alpha=0.5, linewidth=0.75, label=f"Std $\sigma={std:.3f}{dimension}, {std/mean*100:.3f}%$")
            plt.axvline(mean - std, color="black", linestyle="--", alpha=0.5, linewidth=0.75)

            # Print mu and sigma on the plot
            plt.text(mean, np.max(hist), f"$\mu$", ha="center", va="bottom", color="black")
            plt.text(mean + std, np.max(hist), f"$\mu+\sigma$", ha="center", va="bottom", color="black")
            plt.text(mean - std, np.max(hist), f"$\mu-\sigma$", ha="center", va="bottom", color="black")

            # Plot the peaks
            if k == 0:
                coordinate_name = f"{bond_info.atoms[coordinate['i']]['atom']}-{bond_info.atoms[coordinate['j']]['atom']}"
                plt.title(f"Bond length {coordinate_name}: ({mean:.3f} ± {std:.3f})Å\n{len(peaks)} peaks, {normality:.3f}, fixed: {is_fixed}")
            elif k == 1:
                coordinate_name = f"{bond_info.atoms[coordinate['i']]['atom']}-{bond_info.atoms[coordinate['j']]['atom']}-{bond_info.atoms[coordinate['k']]['atom']}"
                plt.title(f"Angle {coordinate_name}: ({mean:.3f} ± {std:.3f})°\n{len(peaks)} peaks, {normality:.3f}, fixed: {is_fixed}")
            elif k == 2:
                coordinate_name = f"{bond_info.atoms[coordinate['i']]['atom']}-{bond_info.atoms[coordinate['j']]['atom']}-{bond_info.atoms[coordinate['k']]['atom']}-{bond_info.atoms[coordinate['l']]['atom']}"
                plt.title(f"Dihedral {coordinate_name}: ({mean:.3f} ± {std:.3f})°\n{len(peaks)} peaks, {normality:.3f}, fixed: {is_fixed}")

            # Axes label
            xLabel = "Bond length" if k == 0 else "Angle" if k == 1 else "Dihedral"
            plt.xlabel(f"{xLabel} [{dimension}]")
            plt.ylabel("Occurances")

            # Plot legend
            plt.legend(fontsize='small')

            name = "bond" if k == 0 else "angle" if k == 1 else "dihedral"
            plt.savefig(f"tmp/plot_{name}_{coordinate_index}.png")
            plt.clf()

            if is_fixed:
                fixed[k].append(coordinate_index)
            else:
                freedom[k].append(coordinate_index)

            coodinates_with_fixed_flag[k].append({"fixed": is_fixed, mean: mean, std: std, **coordinate})

    # Print the results
    print("Bonds:")
    print(f"\tFixed: {len(fixed[0])}")
    print(f"\tFreedom: {len(freedom[0])}")
    print("Angles:")
    print(f"\tFixed: {len(fixed[1])}")
    print(f"\tFreedom: {len(freedom[1])}")
    print("Dihedrals:")
    print(f"\tFixed: {len(fixed[2])}")
    print(f"\tFreedom: {len(freedom[2])}")
    
    print(coodinates_with_fixed_flag[0])
    
    total_dof = len(freedom[0]) + len(freedom[1]) + len(freedom[2])
    total_fixed_dof = len(fixed[0]) + len(fixed[1]) + len(fixed[2])
    print(f"Total degrees of freedom: {total_dof}")
    print(f"{total_dof + total_fixed_dof} -> {total_dof}")