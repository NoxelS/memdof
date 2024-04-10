from typing_extensions import TypedDict

MoleculeInfo = TypedDict(
    "MoleculeInfo",
    {
        "name": str,
        "nrexcl": int,
    },
)

AtomInfo = TypedDict(
    "AtomInfo",
    {
        "nr": int,
        "type": str,
        "resnr": int,
        "residue": str,
        "atom": str,
        "cgnr": int,
        "charge": float,
        "mass": float,
    },
)

BondInfo = TypedDict(
    "BondInfo",
    {
        "i": int,
        "j": int,
        "func": int,
    },
)

PairInfo = TypedDict(
    "PairInfo",
    {
        "i": int,
        "j": int,
        "func": int,
    },
)

AngleInfo = TypedDict(
    "AngleInfo",
    {
        "i": int,
        "j": int,
        "k": int,
        "func": int,
    },
)

DihedralInfo = TypedDict(
    "DihedralInfo",
    {
        "i": int,
        "j": int,
        "k": int,
        "l": int,
        "func": int,
    },
)

ExtendedBondInfo = TypedDict(
    "ExtendedBondInfo",
    {
        "i": int,
        "j": int,
        "func": int,
        "fixed": bool,
        "mean": float,
        "std": float,
    },
)

ExtendedAngleInfo = TypedDict(
    "ExtendedAngleInfo",
    {
        "i": int,
        "j": int,
        "k": int,
        "func": int,
        "fixed": bool,
        "mean": float,
        "std": float,
    },
)

ExtendedDihedralInfo = TypedDict(
    "ExtendedDihedralInfo",
    {
        "i": int,
        "j": int,
        "k": int,
        "l": int,
        "func": int,
        "fixed": bool,
        "mean": float,
        "std": float,
    },
)
