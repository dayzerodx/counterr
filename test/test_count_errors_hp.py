import unittest
from .gen_test_data import gen_mock_hp
from counterr.counters import count_errors_hp
from counterr.util import *

class TestFindHpRegion(unittest.TestCase):
    def test_on_mock_case(self):
        # Generate mock data
        L = int(1e3)
        len_max_hp = 10
        for _ in range(10):
            ref, rd, cig, hp, hist_len_hp_truth = gen_mock_hp(L, len_max_hp, return_hist_len_hp = True)
            alns_scrubbed = [(ref, rd, cig, "", hp)]

            # Run test routine
            hist_len_hp = count_errors_hp(alns_scrubbed, len_min_hp=3, len_max_hp=len_max_hp)

            # Test
            if not (np.allclose(hist_len_hp, hist_len_hp_truth)):
                print(ref)
                print(rd)
                print(cig)
                print(hp)
                # print(hist_len_hp_truth)        
                diff = hist_len_hp_truth - hist_len_hp
                for k in range(4):
                    print(diff[k])
            self.assertTrue(np.allclose(hist_len_hp, hist_len_hp_truth))

        return

if __name__ == "__main__":
    unittest.main()