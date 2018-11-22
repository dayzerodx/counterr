import unittest
from counterr.counters import find_hp_region
from counterr.util import *
from .gen_test_data import gen_mock_hp

class TestFindHpRegion(unittest.TestCase):
    def test_on_mock_case(self):
        # Generate mock data
        L = int(1e5)
        len_max_hp = 20
        ref, rd, cig, hp_truth = gen_mock_hp(L, len_max_hp)

        # Test the function
        hp = find_hp_region(ref, rd, cig, len_min_hp=3)
        for idx_start in range(0, L, 100):
            idx_end = min(L, idx_start+100)
            if hp_truth[idx_start:idx_end] != hp[idx_start:idx_end]:
                print(ref[idx_start:idx_end])
                print(rd[idx_start:idx_end])
                print(cig[idx_start:idx_end])
                print(hp_truth[idx_start:idx_end])
                print(hp[idx_start:idx_end])
                assert False
        self.assertEqual(hp_truth, hp)

        return

if __name__ == "__main__":
    unittest.main()