import os
import unittest
from creeds.matrix_ops import *

class TestJaccardIndex(unittest.TestCase):
	'''
	Test if Jaccard index calculation is correct.
	'''
	def setUp(self):
		## simulate sparse matrices and vectors
		mat_up, vec_up = simulate_matrix(2000, 30000)
		mat_dn, vec_dn = simulate_matrix(2000, 30000)
		self.mat_up = mat_up
		self.mat_dn = mat_dn
		self.vec_up = vec_up
		self.vec_dn = vec_dn

	def test_jaccard(self):
		j = fast_jaccard(self.mat_up, self.vec_up)

		self.assertEquals( j.shape, (2000,) )

		jd = pairwise_distances(self.mat_up.todense(), self.vec_up.reshape(1,-1), 
			metric='jaccard').ravel()
		# Assert allclose with the sklearn implementation
		self.assertTrue(np.allclose(1-j, jd))

		sj = fast_signed_jaccard(self.mat_up, self.mat_dn, 
			self.vec_up, self.vec_dn)

		self.assertEquals(j.shape, sj.shape)

