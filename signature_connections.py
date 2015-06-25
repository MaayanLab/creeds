## 0. get unique signatures across the 3 microtasks using valid_*_entries.pkl
## 1. compute pairwise signed jaccard matrix (issue#4)
## 2. visualize signatures using pack layout or manifold learning (issue#3)
## 3. ...
## created on 6/24/2015


import os, sys
import gzip
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 ## Output Type 3 (Type3) or Type 42 (TrueType)
rcParams['font.sans-serif'] = 'Arial'

import time
import urllib2
import json
import MySQLdb
import cPickle as pickle
from itertools import combinations
from collections import Counter, OrderedDict
from pprint import pprint

from GEOentry import *
sys.path.append('C:\Users\Zichen\Documents\\bitbucket\maayanlab_utils')
from fileIO import file2list, mysqlTable2dict
sys.path.append('C:\Users\Zichen\Documents\GitHub')
from clustergram import clustergram

def list2file(l, fn):
	with open (fn, 'w') as out:
		out.write('\n'.join(map(lambda x :x.encode('utf-8'), l)) + '\n')
	return

# def index_matrix(l, func):
# 	## make index matrix
# 	mat = np.ones((len(l), len(l)))
# 	for i, j in combinations(range(len(l)), 2):
# 		mat[i, j] = func(l[i], l[j])
# 		mat[j, i] = func(l[j], l[i])
# 		if i % 50 == 0 and i != 0:
# 			print i, len(l)
# 	return mat	

def adjMatrx2edgeList(mat, row_labels, col_labels):
	mat_flat = mat.ravel()
	srt_idx = mat_flat.argsort()[::-1]
	edges = []
	for row_label in row_labels:
		for col_label in col_labels:
			edge = row_label + '\t' + col_label
			edges.append(edge)
	edges = np.array(edges)[srt_idx]
	return zip(edges, mat_flat[srt_idx])

def jaccard(l1,l2):
	# case insensitive jaccard
	s1, s2 = set(map(lambda x: x.upper(), l1)), set(map(lambda x:x.upper(), l2))
	up = len(s1 & s2)
	dn = len(s1 | s2)
	if dn == 0: # to prevent 0 division error
		return 0
	else:
		return float(up)/dn


def signed_jaccard(d1, d2):
	## signed jaccard index
	## for dict (overwrite the one for entries)
	j1 = jaccard(d1['up'], d2['up'])
	j2 = jaccard(d1['dn'], d2['dn'])
	j3 = jaccard(d1['dn'], d2['up'])
	j4 = jaccard(d1['up'], d2['dn'])
	return (j1 + j2 - j3 - j4) / 2

def read_gmt_sigs(fn, prefix):
	# read signatures from gmt file
	d = {} # {prefix_id: {'up':[genes], 'dn':[genes]}, ...}
	with open (fn) as f:
		for line in f:
			sl = line.strip().split('\t')
			direction = sl[0].split('|')[1] # up or dn
			uid = sl[1]
			genes = sl[2:]
			prefix_id = prefix + ':' + uid
			if prefix_id not in d:
				d[prefix_id] = {}
			d[prefix_id][direction] = genes
	return d

def get_uniq_sigs():
	# get the uids of signatures that are unique across all microtasks
	all_valid_entry_fns = [
		('dz','output/microtask_dz_jsons/valid_dz_entries.pkl','output/microtask_dz_jsons/crowdsourced_diseases_top600.gmt'),
		('drug','output/microtask_drug_jsons/valid_drug_entries.pkl','output/microtask_drug_jsons/crowdsourced_drugs_top600.gmt'),
		('gene','output/microtask_gene_jsons/valid_gene_entries.pkl','output/microtask_gene_jsons/crowdsourced_single_gene_pert_top600.gmt'),
	] # the order determine priority
	unique_entries = {}
	unique_genesets = {}
	for prefix, valid_entry_fn, gmt_fn in all_valid_entry_fns:
		valid_entries = pickle.load(open(valid_entry_fn, 'rb'))
		d_genesets = read_gmt_sigs(gmt_fn, prefix) 
		for uid, entry in valid_entries.items():
			bools = [valid_e == entry for valid_e in unique_entries.values()]
			prefix_id = prefix + ':' + str(uid) # id to identify entry across microtasks
			if sum(bools) == 0: # all False
				if prefix_id in d_genesets: # some entries failed to convert into geneset possibly due to json encoding
					unique_entries[prefix_id] = entry
					unique_genesets[prefix_id] = d_genesets[prefix_id]

	print 'Number of unique entries across microtasks:',len(unique_entries)
	unique_entries = OrderedDict(sorted(unique_entries.items(), key=lambda t: t[0]))
	unique_genesets = OrderedDict(sorted(unique_genesets.items(), key=lambda t: t[0]))
	return unique_entries, unique_genesets

def pairwise_signed_jaccard(unique_genesets, outfn):
	## compute pairwise signed jaccard and write into file
	with gzip.open (outfn, 'w') as out:
		for i in range(len(unique_genesets)):
			if i % 50 == 0:
				print i, len(unique_genesets)
			prefix_id_i = unique_genesets.keys()[i]
			geneset_i = unique_genesets[prefix_id_i]
			for j in range(i+1, len(unique_genesets)):
				prefix_id_j = unique_genesets.keys()[j]
				geneset_j = unique_genesets[prefix_id_j]
				sj = signed_jaccard(geneset_i, geneset_j)
				out.write('\t'.join(map(str, [prefix_id_i, prefix_id_j, sj] )) + '\n')
	return

unique_entries, unique_genesets = get_uniq_sigs()
# for idx, geneset in unique_genesets.items():
# 	if len(geneset['up']) ==0 or len(geneset['dn']) == 0:
# 		print idx, len(geneset['up']), len(geneset['dn'])

# pairwise_signed_jaccard(unique_genesets, 'signed_jaccard_%s_unique_entries.txt.gz' % len(unique_entries))

## plot embedding for the adjacency matrix
sys.path.append('C:\Users\Zichen\Documents\\bitbucket\\natural_products')
from SparseAdjacencyMat import SparseAdjacencyMat

d_prefix_id_idx = dict(zip(unique_entries.keys(), range(len(unique_entries))))

def read_sam_from_file(fn, d_prefix_id_idx, cutoff=-2):
	## assumes the file is gzipped
	mat = {}
	c = 0
	with gzip.open(fn, 'rb') as f:
		for line in f:
			c += 1
			sl = line.strip().split('\t')
			i, j = d_prefix_id_idx[sl[0]], d_prefix_id_idx[sl[1]]
			score = float(sl[2])
			if score > cutoff:
				mat[i, j] = score
			if c % 2000000 == 0:
				print c
	fn = fn.split('/')[-1]
	print 'finished reading %s' % fn
	return SparseAdjacencyMat(mat, fn)

sam = read_sam_from_file('signed_jaccard_%s_unique_entries.txt.gz' % len(unique_entries), d_prefix_id_idx)
# sam.plot_embedding('truncatedSVD')
sam.plot_embedding('TSNE')

