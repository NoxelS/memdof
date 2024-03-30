from mtypes import *


class BondInfo:

    def __init__(
        self,
        moleculetype: MoleculeInfo,
        atoms: list[AtomInfo],
        bonds: list[BondInfo],
        pairs: list[PairInfo],
        angles: list[AngleInfo],
        dihedrals: list[DihedralInfo],
    ):
        self.moleculetype: MoleculeInfo = moleculetype
        self.atoms: list[AtomInfo] = atoms
        self.bonds: list[BondInfo] = bonds
        self.pairs: list[PairInfo] = pairs
        self.angles: list[AngleInfo] = angles
        self.dihedrals: list[DihedralInfo] = dihedrals


def parse_topology(path: str, ignore_hydrogen: bool = False) -> BondInfo:
    (
        moleculetype_lines,
        atom_lines,
        bond_lines,
        pair_lines,
        angle_lines,
        dihedral_lines,
    ) = ([], [], [], [], [], [])

    with open(path, "r") as f:

        active_directive = None

        for pointer, line in enumerate([l.strip() for l in f.readlines()]):
            # Skip comments and empty lines
            if line.startswith(";") or line.startswith("#") or len(line) == 0:
                continue

            # Switch active directive
            if line.startswith("["):
                active_directive = line.strip("[]").strip().lower()
                continue

            # Parse atoms
            if active_directive == "moleculetype":
                moleculetype_lines.append(line.split())
                continue

            if active_directive == "atoms":
                atom_lines.append(line.split())
                continue

            if active_directive == "bonds":
                bond_lines.append(line.split())
                continue

            if active_directive == "pairs":
                pair_lines.append(line.split())
                continue

            if active_directive == "angles":
                angle_lines.append(line.split())
                continue

            if active_directive == "dihedrals":
                dihedral_lines.append(line.split())
                continue

            # Raise an error if the directive is unknown
            raise ValueError(
                f"Unknown directive '{active_directive}' in line {pointer} of file '{path}'"
            )

    # Parse lines to dictionaries
    moleculetype = [{"name": l[0], "nrexcl": int(l[1])} for l in moleculetype_lines]
    atoms = [
        {
            "nr": int(l[0]) - 1,
            "type": l[1],
            "resnr": int(l[2]),
            "residue": l[3],
            "atom": l[4],
            "cgnr": int(l[5]) - 1,
            "charge": float(l[6]),
            "mass": float(l[7]),
        }
        for l in atom_lines
    ]
    bonds = [
        {"i": int(l[0]) - 1, "j": int(l[1]) - 1, "func": int(l[2])} for l in bond_lines
    ]
    pairs = [
        {"i": int(l[0]) - 1, "j": int(l[1]) - 1, "func": int(l[2])} for l in pair_lines
    ]
    angles = [
        {"i": int(l[0]) - 1, "j": int(l[1]) - 1, "k": int(l[2]) - 1, "func": int(l[3])}
        for l in angle_lines
    ]
    dihedrals = [
        {
            "i": int(l[0]) - 1,
            "j": int(l[1]) - 1,
            "k": int(l[2]) - 1,
            "l": int(l[3]) - 1,
            "func": int(l[4]),
        }
        for l in dihedral_lines
    ]

    moleculetype = moleculetype[0]  # Only one moleculetype is allowed

    if ignore_hydrogen:
        bonds = [
            b
            for b in bonds
            if "H" not in atoms[int(b["i"])]["atom"]
            and "H" not in atoms[int(b["j"])]["atom"]
        ]
        pairs = [
            p
            for p in pairs
            if "H" not in atoms[int(p["i"])]["atom"]
            and "H" not in atoms[int(p["j"])]["atom"]
        ]
        angles = [
            a
            for a in angles
            if "H" not in atoms[int(a["i"])]["atom"]
            and "H" not in atoms[int(a["j"])]["atom"]
            and "H" not in atoms[int(a["k"])]["atom"]
        ]
        dihedrals = [
            d
            for d in dihedrals
            if "H" not in atoms[int(d["i"])]["atom"]
            and "H" not in atoms[int(d["j"])]["atom"]
            and "H" not in atoms[int(d["k"])]["atom"]
            and "H" not in atoms[int(d["l"])]["atom"]
        ]

    return BondInfo(
        moleculetype,
        atoms,
        bonds,
        pairs,
        angles,
        dihedrals,
    )


if __name__ == "__main__":
    pass
