import unittest
from counterr.counters import reconstruct_alignment
from .gen_test_data import gen_rand_seq, MockPysamRead
from counterr.util import *

class TestReconstructAlignment(unittest.TestCase):
	def test_on_mock_read(self):
		# Set numpy random seed
		np.random.seed(42)
		L_contig = int(1e5)

		# Generate reference sequence
		contig = gen_rand_seq(L_contig)

		# Generate mock read
		L_read = int(5e4) 
		idx_start = int(1e4)
		R_ins = 1e-1
		R_del = 1e-1
		R_x = 5e-1
		read = MockPysamRead(contig, L_read, idx_start, R_ins, R_del, R_x)

		# Compute alignment using the function being tested
		ref, rd, cig, Qs = reconstruct_alignment(contig, read, len_max_indel = 20, correct_orientation = True, len_trim_contig_edge=0)

		# Compare to the true answer
		self.assertEqual(ref, read.ref)
		self.assertEqual(rd, read.rd)
		self.assertEqual(cig, read.cig)
		self.assertTrue(np.allclose(Qs, read.Qs))

if __name__ == "__main__":
	unittest.main()