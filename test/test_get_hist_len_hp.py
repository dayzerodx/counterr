import unittest
from counterr.counters import get_hist_len_hp
from counterr.util import *

class TestGetHistLenHp(unittest.TestCase):
    def test_on_mock_contigs(self):
        len_trim_contig_edge = 0
        len_min_hp = 3
        len_max_hp = 10

        # Generate mock contigs and compute the true answer at the same time
        hist_len_hp_truth = np.zeros((N_bases, len_max_hp))
        N_contigs = 100
        L_contig = int(1e3)
        contigs = []
        for _ in range(N_contigs):
            ref = []
            ltr_previous = ""
            idx = 0
            while idx < L_contig:
                ltr = np.random.choice(bases_list)[0]
                if ltr != ltr_previous:
                    ltr_previous = ltr
                    len_hp = np.random.randint(low=len_min_hp, high=len_max_hp, size=1)[0]
                    hist_len_hp_truth[base2int[ltr], len_hp] += 1
                    idx += len_hp
                    ref.append(ltr * len_hp)
            contigs.append(("", "".join(ref)))
        # print(contigs)

        hist_len_hp = get_hist_len_hp(contigs, len_min_hp, len_max_hp, len_trim_contig_edge, verbose=False)
        # Compare to the true answer
        # print(hist_len_hp - hist_len_hp_truth)
        self.assertTrue(np.allclose(hist_len_hp, hist_len_hp_truth))

        return

if __name__ == "__main__":
    unittest.main()