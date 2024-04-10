import pickle
import unittest
from datetime import timedelta
from timeit import default_timer as timer

import numpy as np

from memdof.membranes import calc_IDOF, parse_PTB
from memdof.topologies import parse_topology

TEST_TOPOLOGIE_FILE = "tests/assets/DOPC.itp"


class TestTopologies(unittest.TestCase):

    def test_performance(self):
        start = timer()
        for _ in range(1000):
            bond_info = parse_topology(TEST_TOPOLOGIE_FILE)
        end = timer()
        d1 = timedelta(seconds=(end - start) / 1000)
        self.assertLess(
            d1,
            timedelta(microseconds=4000),
            "Should only take a few milliseconds to parse the topology file.",
        )

    def test_bond_info(self):

        bond_info = parse_topology(TEST_TOPOLOGIE_FILE)
        bond_info_woH = parse_topology(TEST_TOPOLOGIE_FILE, ignore_hydrogen=True)

        self.assertEqual(bond_info.moleculetype["name"], "DOPC")
        self.assertEqual(bond_info_woH.moleculetype["name"], "DOPC")

        self.assertEqual(len(bond_info.atoms), 138)
        self.assertEqual(len(bond_info_woH.atoms), 54)

        self.assertEqual(bond_info.atoms[0]["atom"], "N")
        self.assertEqual(bond_info_woH.atoms[0]["atom"], "N")

    def test_generate_dof():
        # Load the topology information
        topology_info = parse_topology("tests/assets/DOPC.itp", ignore_hydrogen=True)

        # Load the structure of the membrane
        structure, pbc = parse_PTB("tests/assets/3.pdb", "DOPC", quiet=True)

        # Perform the analysis
        extended_topology_info = calc_IDOF(
            structure,
            pbc,
            topology_info,
            quiet=False,
            create_plots=True,
            create_csv=True,
        )

        # Save extended topology info as json
        with open("tmp/extended_toplogy.json", "w") as f:
            f.write(extended_topology_info.json())

        # Save the results as pickle
        with open("tmp/extended_toplogy.pkl", "wb") as f:
            pickle.dump(extended_topology_info, f)


if __name__ == "__main__":
    unittest.main(verbosity=4)
