'''
Implementation of fast Jaccard Index, signed Jaccard Index using sparse matrix operations.
'''
from __future__ import division

import numpy as np
import scipy.sparse as sp


def simulate_matrix(n_rows, n_cols, n_ones_per_row=300, sparse=True):
	'''Simulate a sparse binary matrix and a binary vector. '''
	if not sparse:
		mat = np.zeros((n_rows, n_cols), dtype=np.int8)
	else:
		mat = sp.lil_matrix((n_rows, n_cols), dtype=np.int8)	
	
	for i in xrange(n_rows):
		mat[i, np.random.choice(range(n_cols), n_ones_per_row)] = 1

	if sparse: # convert lil_matrix to csr_matrix for efficient dot product
		mat = mat.tocsr()

	vec = np.zeros(n_cols, dtype=np.int8)
	vec[np.random.choice(range(n_cols), n_ones_per_row)] = 1
	return mat, vec


def fast_jaccard(mat, vec):
	'''Compute row-wise Jaccard index between a matrix and a vector.
	Assumes mat and vec only contain 1's and 0's.
	Return array 
	'''
	assert isinstance(mat, sp.csr_matrix)
	intersections = mat.dot(vec) # number of intersections
	row_sums = mat.getnnz(axis=1).ravel()
	unions = row_sums + vec.sum() - intersections
	return intersections / unions


def fast_signed_jaccard(mat_up, mat_dn, vec_up, vec_dn):
	'''Compute row-wise signed Jaccard between (mat_up, mat_dn) and vec.
	'''
	j1 = fast_jaccard(mat_up, vec_up)
	j2 = fast_jaccard(mat_dn, vec_dn)
	j3 = fast_jaccard(mat_dn, vec_up)
	j4 = fast_jaccard(mat_up, vec_dn)
	return (j1 + j2 - j3 - j4) / 2

