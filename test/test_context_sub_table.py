import unittest
from .gen_test_data import gen_mock_context_sub
from counterr.counters import context_sub_table
from counterr.util import *

class TestContextSub(unittest.TestCase):
    def test_on_mock_case(self):
        # Generate mock data
        L = int(1e3)
        len_context_sub=5
        for _ in range(10):
            ref, rd, cig, sub_matrix_truth = gen_mock_context_sub(L, len_context_sub=len_context_sub)
            alns_scrubbed = [(ref, rd, cig, "", "")]

            # Run test routine
            sub_matrix = context_sub_table(alns_scrubbed, len_context_sub=len_context_sub)

            # Test
            if not np.allclose(sub_matrix, sub_matrix_truth):
                print(ref)
                print(rd)
                print(cig)
                for idx in range(sub_matrix.shape[0]):
                    diff = sub_matrix[idx] - sub_matrix_truth[idx]
                    if diff.sum() != 0:
                        print(idx_to_kmer(len_context_sub, idx))                        
                        print("Truth", sub_matrix_truth[idx])
                        print("Count", sub_matrix[idx])                        
                        print("\n")
            self.assertTrue(np.allclose(sub_matrix, sub_matrix_truth))

        return

if __name__ == "__main__":
    unittest.main()