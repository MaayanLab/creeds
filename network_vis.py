## to make json graph visualizing the signature network
## created on 6/29/2015

## file needed: "signed_jaccard_4066_unique_entries.txt.gz"
## use R to do the hierarchical clustering and 
## use the R package "ape" to retrieve the network edgelist
## from the HC dendrogram


import os, sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn import manifold
import networkx as nx
from networkx.readwrite import json_graph
from sklearn.preprocessing import LabelBinarizer
from sklearn.cluster import AgglomerativeClustering

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 ## Output Type 3 (Type3) or Type 42 (TrueType)
rcParams['font.sans-serif'] = 'Arial'

import json
import MySQLdb
import cPickle as pickle
from itertools import combinations
from collections import Counter
from pprint import pprint
from scipy.stats import pearsonr

HOME = 'C:/Users/Zichen'
sys.path.append(HOME + '/Documents/bitbucket/maayanlab_utils')
from fileIO import read_df, read_gmt, mysqlTable2dict, write_gmt
from plots import COLORS20, COLORS20b





def make_directed_json_graph(gmt_fn, d_id_name, d_id_category, d_category_color, outfn=None):
	# perform HC and make a directed graph and write to json
	# for pack visualization
	d_gmt = read_gmt(gmt_fn)
	d_gmt_filt = {}
	for term, genes in d_gmt.items():
		if len(genes) >= 5:
			d_gmt_filt[term] = genes
	d_gmt = d_gmt_filt

	print 'number of terms:', len(d_gmt)
	umls_ids_kept = d_gmt.keys()
	adj_matrix = jaccard_matrix(d_gmt)

	hc = AgglomerativeClustering(n_clusters=10)
	hc.fit(adj_matrix)

	m = adj_matrix > 0.2
	adj_matrix = adj_matrix * m.astype(int)
	Gu = nx.from_numpy_matrix(adj_matrix) # undirected Graph, to get size

	G = nx.DiGraph()
	for i in range(adj_matrix.shape[0]):
		cluster_label = hc.labels_[i]
		umls_id = umls_ids_kept[i]
		name = d_id_name[umls_id]
		G.add_edge('root', cluster_label)
		G.add_edge(cluster_label, umls_id)
		G.node[umls_id]['size'] = Gu.degree(i)
		G.node[umls_id]['label'] = name

		category = d_id_category[umls_id]
		color = d_category_color[category]
		G.node[umls_id]['color'] = color
	graph_data = json_graph.tree_data(G,root='root')
	json.dump(graph_data, open(outfn, 'wb'))
	return

def make_directed_json_graph_soc(gmt_fn, d_id_name, d_id_category, d_category_color, outfn=None):
	# make directed graph based on SOC - PT
	d_gmt = read_gmt(gmt_fn)
	d_gmt_filt = {}
	for term, genes in d_gmt.items():
		if len(genes) >= 5:
			d_gmt_filt[term] = genes
	d_gmt = d_gmt_filt

	print 'number of terms:', len(d_gmt)
	umls_ids_kept = d_gmt.keys()
	adj_matrix = jaccard_matrix(d_gmt)
	m = adj_matrix > 0.2
	adj_matrix = adj_matrix * m.astype(int)
	Gu = nx.from_numpy_matrix(adj_matrix) # undirected Graph, to get size
	G = nx.DiGraph()
	for i in range(len(umls_ids_kept)):
		umls_id = umls_ids_kept[i]
		name = d_id_name[umls_id]
		category = d_id_category[umls_id]
		color = d_category_color[category]

		G.add_edge('root', category)
		G.add_edge(category, umls_id)

		G.node[umls_id]['size'] = Gu.degree(i)
		G.node[umls_id]['label'] = name
		G.node[umls_id]['color'] = color		
	graph_data = json_graph.tree_data(G,root='root')
	json.dump(graph_data, open(outfn, 'wb'))
	return

