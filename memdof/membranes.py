import csv
import os
import pickle

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.vectors import calc_dihedral
from scipy.signal import find_peaks, savgol_filter
from scipy.stats import normaltest

from memdof.topologies import (ExtendedTopologyInfo, TopologyInfo,
                               parse_topology)

# plt.style.use("seaborn-paper")
# matplotlib.rcParams.update({"font.size": 16})


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


def calc_IDOF(
    structure: Structure,
    pbc: tuple[float, float, float],
    topology_info: TopologyInfo,
    quiet: bool = True,
    max_normality: tuple[float, float, float] = (115, 105, 105),
    max_std: tuple[float, float, float] = (0.1, 1.0, 1.0),
    create_plots: bool = False,
    plots_path: str = "tmp",
    create_csv: bool = False,
    csv_path: str = "tmp",
) -> ExtendedTopologyInfo:
    """
    Calculate the internal degrees of freedom for a given structure.

    Args:
        structure (Structure): The structure object
        pbc (tuple[float, float, float]): The periodic boundary conditions
        topology_info (TopologyInfo): The bond information
        quiet (bool, optional): Suppress warnings. Defaults to True.
        max_normality (tuple[float, float, float], optional): The maximum normality for bonds according to D'Agostino and Pearson's normality test. The tuple contains the maximum normality for bonds, angles and dihedrals. Defaults to (115, 115, 105).
        max_std (tuple[float, float, float], optional): The maximum standard deviation for bonds, angles and dihedrals. Defaults to (0.1 A, 0.5°, 0.5°).
        create_plots (bool, optional): If true creates plots. Defaults to False.
        plots_path (str, optional): The path to save the plots. Defaults to "tmp".
        create_csv (bool, optional): If true creates csv files. Defaults to False.
        csv_path (str, optional): The path to save the csv files. Defaults to "tmp".
    """

    def gaussian(x, mu, sig):
        return (
            1.0
            / (np.sqrt(2.0 * np.pi) * sig)
            * np.exp(-np.power((x - mu) / sig, 2.0) / 2)
        )

    norms = [[] for _ in topology_info.bonds]
    angles = [[] for _ in topology_info.angles]
    dihedrals = [[] for _ in topology_info.dihedrals]

    for chain in structure.get_chains():
        for residue in chain.get_residues():
            atoms = list(residue.get_atoms())

            # Iterate over all bonds
            for coordinate_index, coordinate in enumerate(topology_info.bonds):
                # Get atom names of the bond
                i_name, j_name = (
                    topology_info.atoms[coordinate["i"]]["atom"],
                    topology_info.atoms[coordinate["j"]]["atom"],
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
            for angle_index, angle in enumerate(topology_info.angles):
                # Get atom names of the bond
                i_name, j_name, k_name = (
                    topology_info.atoms[angle["i"]]["atom"],
                    topology_info.atoms[angle["j"]]["atom"],
                    topology_info.atoms[angle["k"]]["atom"],
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
            for dihedral_index, dihedral in enumerate(topology_info.dihedrals):
                # Get atom names of the bond
                i_name, j_name, k_name, l_name = (
                    topology_info.atoms[dihedral["i"]]["atom"],
                    topology_info.atoms[dihedral["j"]]["atom"],
                    topology_info.atoms[dihedral["k"]]["atom"],
                    topology_info.atoms[dihedral["l"]]["atom"],
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
    fixed = [[], [], []]  # Fixed norms, angles and dihedrals
    freedom = [[], [], []]  # Free norms, angles and dihedrals
    coodinates_with_fixed_flag = [[], [], []]

    # Create dir if required
    if create_plots and not os.path.exists(plots_path):
        os.makedirs(plots_path)

    for k, coordinate_list in enumerate(
        [topology_info.bonds, topology_info.angles, topology_info.dihedrals]
    ):
        if not quiet:
            print(
                "Starting analysis of",
                "bonds" if k == 0 else "angles" if k == 1 else "dihedrals",
            )

        for coordinate_index, coordinate in enumerate(coordinate_list):
            # Calculate the mean and standard deviation of the norms
            mean, std = np.mean(coordinates[k][coordinate_index]), np.std(
                coordinates[k][coordinate_index]
            )

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
            is_fixed = normality < max_normality[k] and std < max_std[k]

            if create_plots:
                # Choose different colors for different coordiantes
                color = "black"
                if k == 0:
                    color = "purple"
                elif k == 1:
                    color = "green"
                elif k == 2:
                    color = "orange"

                # Add grid
                plt.grid()

                # Plot the histogram and the peaks bonds of the original data
                plt.hist(
                    coordinates[k][coordinate_index],
                    bins=100,
                    alpha=0.75,
                    color=color,
                    label="Bond length occurances",
                )

                # Plot the smoothed histogram (same x-axis as the original data)
                plt.plot(
                    bin_edges[1:],
                    hist,
                    color=color,
                    label="Smoothed histogram",
                    alpha=0.75,
                )

                # Plot the theoretical normal distribution (continuous)
                x = np.linspace(
                    np.min(coordinates[k][coordinate_index]),
                    np.max(coordinates[k][coordinate_index]),
                    100,
                )
                y = gaussian(x, mean, std)

                # Scale the normal distribution to the histogram
                y = y / np.max(y) * np.max(hist)

                # Choose color for different types of internal coordinates
                color = "b" if is_fixed else "r"
                opacity = 1 if is_fixed else 0.5
                dimension = "Å" if k == 0 else "°"

                # Plot normal distrubtion
                plt.plot(
                    x,
                    y,
                    color=color,
                    linestyle="dotted",
                    label="Normal distribution",
                    alpha=opacity,
                    linewidth=1.5,
                )

                # Plot the mean and standard deviation
                plt.axvline(
                    mean,
                    color="black",
                    linestyle="--",
                    alpha=0.5,
                    linewidth=0.75,
                    label=f"Mean $\mu={mean:.3f}{dimension}$",
                )
                plt.axvline(
                    mean + std,
                    color="black",
                    linestyle="--",
                    alpha=0.5,
                    linewidth=0.75,
                    label=f"Std $\sigma={std:.3f}{dimension}$",
                )
                plt.axvline(
                    mean - std, color="black", linestyle="--", alpha=0.5, linewidth=0.75
                )

                # Print mu and sigma on the plot
                plt.text(
                    mean,
                    np.max(hist),
                    f"$\mu$",
                    ha="center",
                    va="bottom",
                    color="black",
                )
                plt.text(
                    mean + std,
                    np.max(hist),
                    f"$\mu+\sigma$",
                    ha="center",
                    va="bottom",
                    color="black",
                )
                plt.text(
                    mean - std,
                    np.max(hist),
                    f"$\mu-\sigma$",
                    ha="center",
                    va="bottom",
                    color="black",
                )

                # Plot the peaks
                if k == 0:
                    coordinate_name = f"{topology_info.atoms[coordinate['i']]['atom']}-{topology_info.atoms[coordinate['j']]['atom']}"
                    plt.title(
                        f"Bond length {coordinate_name}\nNormality: {normality:.3f}"
                    )
                elif k == 1:
                    coordinate_name = f"{topology_info.atoms[coordinate['i']]['atom']}-{topology_info.atoms[coordinate['j']]['atom']}-{topology_info.atoms[coordinate['k']]['atom']}"
                    plt.title(f"Angle {coordinate_name}\nNormality: {normality:.3f}")
                elif k == 2:
                    coordinate_name = f"{topology_info.atoms[coordinate['i']]['atom']}-{topology_info.atoms[coordinate['j']]['atom']}-{topology_info.atoms[coordinate['k']]['atom']}-{topology_info.atoms[coordinate['l']]['atom']}"
                    plt.title(f"Dihedral {coordinate_name}\nNormality: {normality:.3f}")

                # Axes label
                xLabel = "Bond length" if k == 0 else "Angle" if k == 1 else "Dihedral"
                plt.xlabel(f"{xLabel} [{dimension}]")
                plt.ylabel("Occurances")

                # Plot legend
                plt.legend(fontsize="small")

                name = "bond" if k == 0 else "angle" if k == 1 else "dihedral"

                if not os.path.exists(os.path.join(plots_path, "fixed")):
                    os.makedirs(os.path.join(plots_path, "fixed"))

                if not os.path.exists(os.path.join(plots_path, "free")):
                    os.makedirs(os.path.join(plots_path, "free"))

                path = (
                    os.path.join(plots_path, "fixed")
                    if is_fixed
                    else os.path.join(plots_path, "free")
                )
                plt.savefig(os.path.join(path, f"{name}_{coordinate_index}.png"))
                plt.clf()

            if is_fixed:
                fixed[k].append(coordinate_index)
            else:
                freedom[k].append(coordinate_index)

            coodinates_with_fixed_flag[k].append(
                {"fixed": is_fixed, "mean": mean, "std": std, **coordinate}
            )

    if not quiet:
        total_dof = len(freedom[0]) + len(freedom[1]) + len(freedom[2])
        total_fixed_dof = len(fixed[0]) + len(fixed[1]) + len(fixed[2])

        # Print the results
        print("Bonds:")
        print(
            f" - Fixed: {len(fixed[0])} ({len(fixed[0])/len(topology_info.bonds)*100:.2f}%)"
        )
        print(
            f" - Free:  {len(freedom[0])} ({len(freedom[0])/len(topology_info.bonds)*100:.2f}%)"
        )
        print("Angles:")
        print(
            f" - Fixed: {len(fixed[1])} ({len(fixed[1])/len(topology_info.angles)*100:.2f}%)"
        )
        print(
            f" - Free:  {len(freedom[1])} ({len(freedom[1])/len(topology_info.angles)*100:.2f}%)"
        )
        print("Dihedrals:")
        print(
            f" - Fixed: {len(fixed[2])} ({len(fixed[2])/len(topology_info.dihedrals)*100:.2f}%)"
        )
        print(
            f" - Free:  {len(freedom[2])} ({len(freedom[2])/len(topology_info.dihedrals)*100:.2f}%)"
        )
        print(f"Total free degrees of freedom:  {total_dof}")
        print(f"Total fixed degrees of freedom: {total_fixed_dof}")
        print(
            f"{total_dof + total_fixed_dof} -> {total_dof} (decreased by {(1-total_dof/(total_dof + total_fixed_dof))*100:.2f}%)"
        )

    if create_csv:
        # Create dir if required
        if not os.path.exists(csv_path):
            os.makedirs(csv_path)

        # Generate paths
        bond_path = os.path.join(csv_path, "bonds.csv")
        angle_path = os.path.join(csv_path, "angles.csv")
        dihedral_path = os.path.join(csv_path, "dihedrals.csv")

        # Make csv for fixed and free coordinates for each type of internal coordinate
        with open(bond_path, "w") as f:
            writer = csv.writer(f, lineterminator="\n", delimiter=";")
            writer.writerow(["i", "j", "fixed", "mean", "std"])
            for i, bond in enumerate(coodinates_with_fixed_flag[0]):
                writer.writerow(
                    [
                        bond["i"],
                        bond["j"],
                        bond["fixed"],
                        f"{bond['mean']:.3f}",
                        f"{bond['std']:.3f}",
                    ]
                )

        with open(angle_path, "w") as f:
            writer = csv.writer(f, lineterminator="\n", delimiter=";")
            writer.writerow(["i", "j", "k", "fixed", "mean", "std"])
            for i, angle in enumerate(coodinates_with_fixed_flag[1]):
                writer.writerow(
                    [
                        angle["i"],
                        angle["j"],
                        angle["k"],
                        angle["fixed"],
                        f"{angle['mean']:.3f}",
                        f"{angle['std']:.3f}",
                    ]
                )

        with open(dihedral_path, "w") as f:
            writer = csv.writer(f, lineterminator="\n", delimiter=";")
            writer.writerow(["i", "j", "k", "l", "fixed", "mean", "std"])
            for i, dihedral in enumerate(coodinates_with_fixed_flag[2]):
                writer.writerow(
                    [
                        dihedral["i"],
                        dihedral["j"],
                        dihedral["k"],
                        dihedral["l"],
                        dihedral["fixed"],
                        f"{dihedral['mean']:.3f}",
                        f"{dihedral['std']:.3f}",
                    ]
                )

    # Build the extended topology info
    extended_topology_info = ExtendedTopologyInfo(
        moleculetype=topology_info.moleculetype,
        atoms=topology_info.atoms,
        pairs=topology_info.pairs,
        bonds=coodinates_with_fixed_flag[0],
        angles=coodinates_with_fixed_flag[1],
        dihedrals=coodinates_with_fixed_flag[2],
    )

    if create_plots:
        """
        Create an image with three plots for bonds, angles and dihedrals.
        Plot all the histograms and the peaks of the original data into the plot.
        """

        fig, ax = plt.subplots(3, 1, figsize=(13, 10))
        ax[0].set_title("Bonds")
        ax[1].set_title("Angles")
        ax[2].set_title("Dihedrals")

        x_ranges = [(99, 0), (180, 0), (180, 0)]

        for k, coordinate_list in enumerate(
            [topology_info.bonds, topology_info.angles, topology_info.dihedrals]
        ):
            for coordinate_index, coordinate in enumerate(coordinate_list):
                # Choose color for different types of internal coordinates
                color = "purple" if coordinate_index in fixed[k] else "red"
                dimension = "Å" if k == 0 else "°"
                zorder = 1 if coordinate_index in fixed[k] else 0
                label = "Fixed" if coordinate_index in fixed[k] else "Free"
                min_alpha = 0.25 if k == 0 else 0.85
                max_alpha = 0.5 if k == 0 else 0.25
                alpha = min_alpha if coordinate_index in fixed[k] else max_alpha

                # Calculate the range of the x-axis
                max_range = np.max(coordinates[k][coordinate_index])
                min_range = np.min(coordinates[k][coordinate_index])

                # Update the range if required
                if min_range < x_ranges[k][0]:
                    x_ranges[k] = (min_range, x_ranges[k][1])
                if max_range > x_ranges[k][1]:
                    x_ranges[k] = (x_ranges[k][0], max_range)

                # Plot the histogram and the peaks bonds of the original data
                ax[k].hist(
                    coordinates[k][coordinate_index],
                    bins=150,
                    alpha=alpha,
                    color=color,
                    zorder=zorder,
                )

                # Add empty plot to add legend
                ax[k].plot([], [], color="purple", alpha=1, label="Fixed")
                ax[k].plot([], [], color="red", alpha=1, label="Free")

                # Axes label
                xLabel = "Bond length" if k == 0 else "Angle" if k == 1 else "Dihedral"
                ax[k].set_xlabel(f"{xLabel} [{dimension}]")
                ax[k].set_ylabel("Occurances")

        # Explain red = free, purple = fixed on the plot as legend
        for i in range(3):
            ax[i].legend(["Fixed", "Free"])

        # Add more ticks on x-axis
        ax[0].set_xticks(
            np.linspace(
                np.floor(x_ranges[0][0] * 10) * 0.1,
                np.ceil(x_ranges[0][1] * 10) * 0.1,
                num=10,
            )
        )
        ax[1].set_xticks(
            np.linspace(
                np.floor(x_ranges[1][0]),
                np.ceil(x_ranges[1][1]),
                num=int((np.ceil(x_ranges[1][1]) - np.floor(x_ranges[1][0])) / 5),
            )
        )
        ax[2].set_xticks(
            np.linspace(
                np.floor(x_ranges[2][0]),
                np.ceil(x_ranges[2][1]),
                num=int((np.ceil(x_ranges[2][1]) - np.floor(x_ranges[2][0])) / 6),
            )
        )

        # Tight layout
        plt.tight_layout()
        plt.savefig(os.path.join(plots_path, "plot_all.png"))
        plt.clf()

    return extended_topology_info
