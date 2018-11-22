import unittest
from .gen_test_data import gen_mock_errors
from counterr.counters import count_errors
from counterr.util import *

class TestCountErrors(unittest.TestCase):
    def test_on_mock_case(self):
        # Generate mock data
        L = int(1e2)
        len_max_indel = 20
        for _ in range(100):
            ref, rd, cig, sub_matrix_truth, hist_len_del_truth, hist_len_ins_truth = gen_mock_errors(L, len_max_indel)
            alns_scrubbed = [(ref, rd, cig, "", "")]

            # Run test routine
            sub_matrix, hist_len_del, hist_len_ins = count_errors(alns_scrubbed, len_min_hp=3, exclude_hp = False, len_max_indel=len_max_indel)

            # Test
            if not (np.allclose(sub_matrix, sub_matrix_truth)
                and np.allclose(hist_len_del, hist_len_del_truth)
                and np.allclose(hist_len_ins, hist_len_ins_truth)):
                print(ref)
                print(rd)
                print(cig)
                print(sub_matrix - sub_matrix_truth)
                print(hist_len_del - hist_len_del_truth)
                print(hist_len_ins - hist_len_ins_truth)
            self.assertTrue(np.allclose(sub_matrix, sub_matrix_truth))
            self.assertTrue(np.allclose(hist_len_del, hist_len_del_truth))
            self.assertTrue(np.allclose(hist_len_ins, hist_len_ins_truth))

        return

if __name__ == "__main__":
    unittest.main()