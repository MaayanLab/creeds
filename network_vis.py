## to make json graph visualizing the signature network
## created on 6/29/2015

## file needed: "signed_jaccard_4066_unique_entries.txt.gz"
## use R to do the hierarchical clustering and 
## use the R package "ape" to retrieve the network edgelist
## from the HC dendrogram


import os, sys
import gzip
import numpy as np
import matplotlib.pyplot as plt
from sklearn import manifold
import networkx as nx
from networkx.readwrite import json_graph


import rpy2.robjects as ro
ro.r('''
	source('HC2graph.R')
''')
hc2network = ro.globalenv['hc2network'] ## convert R function to python function

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 ## Output Type 3 (Type3) or Type 42 (TrueType)
rcParams['font.sans-serif'] = 'Arial'

import json
import MySQLdb
import cPickle as pickle
from itertools import combinations
from collections import Counter, OrderedDict
from pprint import pprint


HOME = 'C:/Users/Zichen'
sys.path.append(HOME + '/Documents/bitbucket/maayanlab_utils')
from fileIO import read_df, read_gmt, mysqlTable2dict, write_gmt
from plots import COLORS10,COLORS20, COLORS20b

#### functions copied from 'signature_connections.py'
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

unique_entries, _ = get_uniq_sigs()

sys.path.append('C:\Users\Zichen\Documents\\bitbucket\\natural_products')
from SparseAdjacencyMat import SparseAdjacencyMat

d_prefix_id_idx = dict(zip(unique_entries.keys(), range(len(unique_entries))))

def read_sam_from_file(fn, d_prefix_id_idx, cutoff=-2):
	## assumes the file is gzipped
	## return a sam whose index are arbituray id encoded in d_prefix_id_idx
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

#### end of functions copied from 'signature_connections.py'


## get some meta-data for the signatures
d_prefix_id_label = {}
d_prefix_id_color = {}
d_prefix_id_size = {} # number of samples
d_dz_name = mysqlTable2dict('maaya0_crowdsourcing', 'geo2enrichr_dz', -1, 3)
d_gene_name = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_genes', 0, 1)
d_drug_name = mysqlTable2dict('maaya0_crowdsourcing', 'geo2enrichr_drug', -1, 3)

for prefix_id in d_prefix_id_idx:
	uid = int(prefix_id.split(':')[1])
	d_prefix_id_size[prefix_id] = len(unique_entries[prefix_id])
	if prefix_id.startswith('dz'):
		d_prefix_id_label[prefix_id] = d_dz_name[uid].lower()
		d_prefix_id_color[prefix_id] = COLORS10[0][1:]
	elif prefix_id.startswith('drug'):
		d_prefix_id_label[prefix_id] = d_drug_name[uid].lower()
		d_prefix_id_color[prefix_id] = COLORS10[1][1:]
	else:
		d_prefix_id_label[prefix_id] = d_gene_name[uid]
		d_prefix_id_color[prefix_id] = COLORS10[2][1:]


def np2matrix(mat):
	## convert a numpy 2D array to R matrix
	from rpy2.robjects import numpy2ri 
	numpy2ri.activate()
	nr, nc = mat.shape
	return ro.r.matrix(mat, nrow=nr, ncol=nc)

def rmatrix2np(rmat):
	## convert an R matrix to numpy array
	return np.array(rmat) - 1 # R is 1 indexed, python is 0 indexed

def make_directed_json_graph(fn, depth=4, outfn=None):
	sam = read_sam_from_file(fn, d_prefix_id_idx)
	adj_matrix = sam.to_csr_matrix().toarray()
	del sam
	adj_matrix = np2matrix(adj_matrix)
	edgelist = hc2network(adj_matrix, 'cosine', 'average')
	edgelist = rmatrix2np(edgelist) # convert to numpy array
	print 'finished generating edgelist'

	prefix_ids = unique_entries.keys() # ordered unique prefix_ids
	all_idx = set(d_prefix_id_idx.values())

	## convert the edgelist to DiGraph
	G = nx.DiGraph()
	for n1, n2 in edgelist:
		G.add_edge(n1, n2)
	del edgelist
	print 'edgelist converted to G'	
	## find root:
	for n,d in G.in_degree().items():
		if d == 0:
			root = n
			break
	print 'root found:', root
	## trim internal branches to shorten the depth of the dendrogram
	G_trimed = nx.DiGraph()
	for leaf in all_idx: # all the leafs 
		shortest_path = nx.shortest_path(G, root, leaf)
		leaf_prefix_id = prefix_ids[leaf] # prefix id of leaf
		leaf_label = d_prefix_id_label[leaf_prefix_id]
		leaf_color = d_prefix_id_color[leaf_prefix_id]
		leaf_size  = d_prefix_id_size[leaf_prefix_id]

		if len(shortest_path) > depth:
			shortest_path = shortest_path[0:depth-1] + [leaf_prefix_id] # remove intermediate nodes
		else: # want the id of leaf node to be prefix_id
			shortest_path[-1] = leaf_prefix_id

		G_trimed.add_path(shortest_path)
		G_trimed.node[leaf_prefix_id]['label'] = leaf_label
		G_trimed.node[leaf_prefix_id]['color'] = leaf_color
		G_trimed.node[leaf_prefix_id]['size'] = leaf_size
		
	del G
	print G_trimed.number_of_edges(), G_trimed.number_of_nodes()
	graph_data = json_graph.tree_data(G_trimed,root=root)
	json.dump(graph_data, open(outfn, 'wb'))
	return

def make_directed_json_graph_cat(outfn=None):
	# make directed graph based on category
	prefix_ids = unique_entries.keys() # ordered unique prefix_ids
	G = nx.DiGraph()
	# for i in range(len(umls_ids_kept)):
	for prefix_id in prefix_ids:
		label = d_prefix_id_label[prefix_id]
		category = prefix_id.split(':')[0]
		color = d_prefix_id_color[prefix_id]
		size = d_prefix_id_size[prefix_id]
		G.add_path(['root', category, prefix_id])

		G.node[prefix_id]['label'] = label
		G.node[prefix_id]['color'] = color
		G.node[prefix_id]['size'] = size
	graph_data = json_graph.tree_data(G,root='root')
	json.dump(graph_data, open(outfn, 'wb'))
	return

### testing
# mat = np.arange(100).reshape(10,10)
# rmat = np2matrix(mat)
# edgelist = hc2network(rmat, 'cosine', 'average')
# # print edgelist
# print np.array(edgelist)

### making json graphs
make_directed_json_graph('signed_jaccard_4066_unique_entries.txt.gz', 
	depth=6, outfn='signature_digraph_depth6.json')

# make_directed_json_graph_cat(outfn='signature_digraph_cat.json')

# print d_prefix_id_idx['gene:2165']
