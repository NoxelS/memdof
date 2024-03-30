import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.vectors import calc_dihedral
from mtypes import *
from scipy.signal import find_peaks, savgol_filter


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

    CUTOFF_STD = 0.01  # Å

    import matplotlib.pyplot as plt
    from topologies import parse_topology

    bond_info = parse_topology("tests/assets/DOPC.itp", ignore_hydrogen=True)

    norms = [[] for _ in bond_info.bonds]
    angles = [[] for _ in bond_info.angles]
    dihedrals = [[] for _ in bond_info.dihedrals]

    for chain in structure.get_chains():
        for residue in chain.get_residues():
            atoms = list(residue.get_atoms())

            # Iterate over all bonds
            for bond_index, bond in enumerate(bond_info.bonds):
                # Get atom names of the bond
                i_name, j_name = (
                    bond_info.atoms[bond["i"]]["atom"],
                    bond_info.atoms[bond["j"]]["atom"],
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
                norms[bond_index].append(np.linalg.norm(vec))

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

    bonds_freedom = []
    bonds_fixed = []

    angles_freedom = []
    angles_fixed = []

    dihedrals_freedom = []
    dihedrals_fixed = []

    for bond_index, bond in enumerate(bond_info.bonds):
        # Calculate the mean and standard deviation of the norms
        mean, std = np.mean(norms[bond_index]), np.std(norms[bond_index])

        # Bin the norms in 100 bins
        hist, bin_edges = np.histogram(norms[bond_index], bins=100)

        # Find value where the histogram is at its peak
        peak0 = np.argmax(hist)

        # Remove all norms that are less than 10% of the peaks occurrence
        hist[hist < 0.05 * hist[peak0]] = 0

        # Apply a Savitzky-Golay filter to the histogram
        hist = savgol_filter(hist, 21, 3)

        # Remove negative values
        hist[hist < 0] = 0

        # Find peaks for histogram
        peaks, dict = find_peaks(
            hist,
            threshold=[0, 99999],
            distance=20,
            width=10,
            height=[np.max(hist) * 0.05, 99999],
        )

        # Plot the histogram and the peaks bonds of the original data
        plt.hist(norms[bond_index], bins=100, alpha=0.5, color="b")

        # Plot the smoothed histogram (same x-axis as the original data)
        plt.plot(bin_edges[1:], hist, color="purple")

        # Plot the mean and standard deviation
        plt.axvline(mean, color="r", linestyle="-", alpha=0.5)
        plt.axvline(mean + std, color="r", linestyle="-", alpha=0.5)
        plt.axvline(mean - std, color="r", linestyle="-", alpha=0.5)

        # Plot the peaks
        plt.plot(bin_edges[peaks], hist[peaks], "x", color="r")
        plt.title(f"{bond_index}: {mean:.3f} ± {std:.3f} Å ({len(peaks)} peaks)")

        plt.savefig(f"tmp/plot_{bond_index}.png")
        plt.clf()

        # Print the results
        # print(f"{bond_index}: {mean:.3f} ± {std:.3f} Å ({len(peaks)} peaks)")

        # TODO: Use more criteria

        if len(peaks) == 1:
            bonds_fixed.append(bond_index)
        else:
            bonds_freedom.append(bond_index)

    for angle_index, angle in enumerate(bond_info.angles):
        # Calculate the mean and standard deviation of the angles
        mean, std = np.mean(angles[angle_index]), np.std(angles[angle_index])

        # Bin the norms in 100 bins
        hist, bin_edges = np.histogram(angles[angle_index], bins=100)

        # Find value where the histogram is at its peak
        peak0 = np.argmax(hist)

        # Remove all norms that are less than 10% of the peaks occurrence
        hist[hist < 0.05 * hist[peak0]] = 0

        # Apply a Savitzky-Golay filter to the histogram
        hist = savgol_filter(hist, 21, 3)

        # Remove negative values
        hist[hist < 0] = 0

        # Find peaks for histogram
        peaks, dict = find_peaks(
            hist,
            threshold=[0, 99999],
            distance=20,
            width=10,
            height=[np.max(hist) * 0.05, 99999],
        )

        # Plot the histogram and the peaks bonds of the original data
        plt.hist(angles[angle_index], bins=100, alpha=0.5, color="b")

        # Plot the smoothed histogram (same x-axis as the original data)
        plt.plot(bin_edges[1:], hist, color="purple")

        # Plot the mean and standard deviation
        plt.axvline(mean, color="r", linestyle="-", alpha=0.5)
        plt.axvline(mean + std, color="r", linestyle="-", alpha=0.5)
        plt.axvline(mean - std, color="r", linestyle="-", alpha=0.5)

        # Plot the peaks
        plt.plot(bin_edges[peaks], hist[peaks], "x", color="r")
        plt.title(f"{angle_index}: {mean:.3f} ± {std:.3f} Å ({len(peaks)} peaks)")

        plt.savefig(f"tmp/plot_angle_{angle_index}.png")
        plt.clf()

        # Print the results
        # print(f"{bond_index}: {mean:.3f} ± {std:.3f} Å ({len(peaks)} peaks)")

        # TODO: Use more criteria

        if len(peaks) == 1:
            angles_fixed.append(bond_index)
        else:
            angles_freedom.append(bond_index)

    for dihedral_index, dihedral in enumerate(bond_info.dihedrals):
        # Calculate the mean and standard deviation of the dihedrals
        mean, std = np.mean(dihedrals[dihedral_index]), np.std(
            dihedrals[dihedral_index]
        )

        # Bin the norms in 100 bins
        hist, bin_edges = np.histogram(dihedrals[dihedral_index], bins=100)

        # Find value where the histogram is at its peak
        peak0 = np.argmax(hist)

        # Remove all norms that are less than 10% of the peaks occurrence
        hist[hist < 0.05 * hist[peak0]] = 0

        # Apply a Savitzky-Golay filter to the histogram
        hist = savgol_filter(hist, 21, 3)

        # Remove negative values
        hist[hist < 0] = 0

        # Find peaks for histogram
        peaks, dict = find_peaks(
            hist,
            threshold=[0, 99999],
            distance=20,
            width=10,
            height=[np.max(hist) * 0.05, 99999],
        )

        # Plot the histogram and the peaks bonds of the original data
        plt.hist(dihedrals[dihedral_index], bins=100, alpha=0.5, color="b")

        # Plot the smoothed histogram (same x-axis as the original data)
        plt.plot(bin_edges[1:], hist, color="purple")

        # Plot the mean and standard deviation
        plt.axvline(mean, color="r", linestyle="-", alpha=0.5)
        plt.axvline(mean + std, color="r", linestyle="-", alpha=0.5)
        plt.axvline(mean - std, color="r", linestyle="-", alpha=0.5)

        # Plot the peaks
        plt.plot(bin_edges[peaks], hist[peaks], "x", color="r")
        plt.title(f"{dihedral_index}: {mean:.3f} ± {std:.3f} Å ({len(peaks)} peaks)")

        plt.savefig(f"tmp/plot_dihedral_{dihedral_index}.png")
        plt.clf()

        # Print the results
        # print(f"{bond_index}: {mean:.3f} ± {std:.3f} Å ({len(peaks)} peaks)")

        # TODO: Use more criteria

        if len(peaks) == 1:
            dihedrals_fixed.append(bond_index)
        else:
            dihedrals_freedom.append(bond_index)

    print(f"Fixed bonds: {len(bonds_fixed)}")
    print(f"Free bonds: {len(bonds_freedom)}")
    print(f"Fixed angles: {len(angles_fixed)}")
    print(f"Free angles: {len(angles_freedom)}")
    print(f"Fixed dihedrals: {len(dihedrals_fixed)}")
    print(f"Free dihedrals: {len(dihedrals_freedom)}")
    
    total_dof = len(bonds_freedom) + len(angles_freedom) + len(dihedrals_freedom)
    total_fixed_dof = len(bonds_fixed) + len(angles_fixed) + len(dihedrals_fixed)
    print(f"Total degrees of freedom: {total_dof}")
    print(f"{total_dof + total_fixed_dof} -> {total_dof}")
    
