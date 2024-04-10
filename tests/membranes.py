import pickle
import unittest
from datetime import timedelta
from timeit import default_timer as timer

import numpy as np

from memdof import calc_IDOF, parse_PTB, parse_topology

TEST_TOPOLOGIE_FILE = "tests/assets/DOPC.itp"


class TestMembranes(unittest.TestCase):

    def test_generate_dof(self):
        # Load the topology information
        topology_info = parse_topology("tests/assets/DOPC.itp", ignore_hydrogen=True)

        # Load the structure of the membrane
        structure, pbc = parse_PTB("tests/assets/1.pdb", "DOPC", quiet=True)

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
