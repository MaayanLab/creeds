## 0. get unique signatures across the 3 microtasks using valid_*_entries.pkl
## 1. compute pairwise signed jaccard matrix (issue#4)
## 2. visualize signatures using pack layout or manifold learning (issue#3)
## 3. ...
## created on 6/24/2015


import os, sys
import gzip
import numpy as np
import networkx as nx
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
from sklearn.metrics import auc
from statsmodels.sandbox.stats.multicomp import multipletests

from GEOentry import *
sys.path.append('C:\Users\Zichen\Documents\\bitbucket\maayanlab_utils')
from fileIO import file2list, mysqlTable2dict
from plots import enlarge_tick_fontsize

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

from pymongo import MongoClient
def get_limma_sigs_from_db(unique_entries, correction_method='fdr_bh'):
	# correction_method should be in ['fdr_bh', 'bonferroni']
	unique_genesets = {}
	cutoff = 0.05
	# client = MongoClient('mongodb://127.0.0.1:27017/')
	client = MongoClient('mongodb://146.203.54.131:27017/')
	db = client['microtask_signatures']
	COLL = db['signatures']
	projection = {'_id':False, 'limma':True, 'fold_changes':True}
	for uid in unique_entries:
		doc = COLL.find_one({'id':uid}, projection)
		if doc is not None:
			if 'limma' in doc and 'fold_changes' in doc:
				d_gene_fc = dict(zip(doc['fold_changes']['genes'], doc['fold_changes']['vals'])) # raw fold change
				pvals = np.array(doc['limma']['vals'])
				genes = np.array(doc['limma']['genes'])
				_, qvals, _, _ = multipletests(pvals, method=correction_method)
				unique_genesets[uid] = {'up':[], 'dn':[]}
				for gene in genes[qvals < cutoff].tolist():
					if gene in d_gene_fc:
						if d_gene_fc[gene] > 1:
							unique_genesets[uid]['up'].append(gene)
						else:
							unique_genesets[uid]['dn'].append(gene)
	return unique_genesets


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

from scipy.spatial.distance import cosine
def chdir_cosine(chdir1, chdir2):
	chdir1 = dict(zip(map(lambda x:x.upper(), chdir1['genes']), chdir1['vals']))
	chdir2 = dict(zip(map(lambda x:x.upper(), chdir2['genes']), chdir2['vals']))

	shared_genes = set(chdir1.keys()) & set(chdir2.keys())
	# shared_genes = list(shared_genes)
	vals1 = [ chdir1[gene] for gene in shared_genes ]
	vals2 = [ chdir2[gene] for gene in shared_genes ]
	return 1 - cosine(vals1, vals2) # cosine similarity


def pairwase_cosine(unique_entries, outfn):
	## calculate cosine distance between chdir signatures in the db
	## by using the shared genes
	client = MongoClient('mongodb://146.203.54.131:27017/')
	db = client['microtask_signatures']
	COLL = db['signatures']
	projection = {'_id':False, 'chdir':True}
	all_uids = unique_entries.keys()
	with gzip.open (outfn, 'w') as out:
		for i, uid_i in enumerate(all_uids):
			doc_i = COLL.find_one({'id':uid_i}, projection)
			if doc_i is not None:
				chdir_i = doc_i['chdir']
				for j in range(i+1, len(all_uids)):
					uid_j = all_uids[j]
					doc_j = COLL.find_one({'id':uid_j}, projection)
					if doc_j is not None:
						chdir_j = doc_j['chdir']
						score = chdir_cosine(chdir_i, chdir_j)
						out.write('\t'.join(map(str, [uid_i, uid_j, score] )) + '\n')
	return	


unique_entries, unique_genesets = get_uniq_sigs()
# for idx, geneset in unique_genesets.items():
# 	if len(geneset['up']) ==0 or len(geneset['dn']) == 0:
# 		print idx, len(geneset['up']), len(geneset['dn'])

# pairwise_signed_jaccard(unique_genesets, 'signed_jaccard_%s_gene_unique_entries.txt.gz' % len(unique_entries))

# unique_genesets = get_limma_sigs_from_db(unique_entries, correction_method='bonferroni')
# pairwase_cosine(unique_entries)
# json.dump(unique_genesets, open('limma_genesets_bonferroni_0.05.json', 'wb'))
# pairwise_signed_jaccard(unique_genesets, 'signed_jaccard_%s_gene_unique_entries_limma_bonferroni_0.05.txt.gz' % len(unique_entries))

pairwase_cosine(unique_entries, 'cosine_similarity_%s_gene_unique_entries.txt.gz' % len(unique_entries))

'''
d_gds_gse = mysqlTable2dict('maaya0_crowdsourcing', 'geo2enrichr_dz', 0, 1)
## load dz id
# global d_uid_doid, d_uid_umls, d_uid_gse
# d_uid_doid = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_dzs', 0, 1)
# d_uid_umls = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_dzs', 0, 2)
# d_uid_geoid = mysqlTable2dict('maaya0_crowdsourcing', 'geo2enrichr_dz', -1, 0)

## load cleaned gene symbols
# global d_uid_hs, d_uid_mm, d_uid_gse
# d_uid_hs = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_genes', 0, 1)
# d_uid_mm = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_genes', 0, 2)
# d_uid_geoid = mysqlTable2dict('maaya0_crowdsourcing', 'geo2enrichr', -1, 0)

## load cleaned drugs
global d_uid_gse
# d_uid_hs = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_genes', 0, 1)
# d_uid_mm = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_genes', 0, 2)
d_uid_geoid = mysqlTable2dict('maaya0_crowdsourcing', 'geo2enrichr_drug', -1, 0)

d_uid_gse = {} # convert gds to gse
for uid, geoid in d_uid_geoid.items():
	if geoid in d_gds_gse: gse = d_gds_gse[geoid]
	else: gse = geoid
	d_uid_gse[uid] = gse
'''

def read_pairwise_from_file(fn, absolute=False):
	## read signed jaccard file
	## return sorted_list
	mat = {}
	id_pairs = []
	scores = []
	with gzip.open(fn, 'rb') as f:
		for line in f:
			sl = line.strip().split('\t')
			uid_i, uid_j = map(lambda x: int(x.split(':')[1]), sl[:2])
			score = float(sl[2])
			id_pairs.append((uid_i, uid_j))
			scores.append(score)
	# sort:
	id_pairs = np.array(id_pairs)
	scores = np.array(scores)
	if absolute:
		srt_idx = np.argsort(abs(scores))[::-1] # large to small
	else:
		srt_idx = np.argsort(scores)[::-1] # large to small
	return id_pairs[srt_idx]

def _plot_roc_same_dz(id_pairs):
	## compute TPR, FPR, auroc for same dz
	fp_list = []
	tp_list = []
	fp = 0
	tp = 0
	for uid_i, uid_j in id_pairs:
		gse_i, gse_j = d_uid_gse[uid_i], d_uid_gse[uid_j]
		if gse_i != gse_j: # exclude signatures from the same GEO study
			doid_i, doid_j = d_uid_doid[uid_i], d_uid_doid[uid_j]
			umls_i, umls_j = d_uid_umls[uid_i], d_uid_umls[uid_j]
			if doid_i == doid_j and None not in [doid_i, doid_j]:
				tp += 1
			elif umls_i == umls_j and None not in [umls_i, umls_j]:
				tp += 1
			else:
				fp += 1
			tp_list.append(tp)
			fp_list.append(fp)
	max_TP = float(tp_list[-1])
	max_FP = float(fp_list[-1])
	TPRs = np.array(tp_list)/max_TP
	FPRs = np.array(fp_list)/max_FP
	auroc = auc(FPRs, TPRs)
	return auroc, FPRs, TPRs

def load_do_network():
	# load disease ontology
	# G is a DiGraph from parent to children
	import obo_parser
	DO_FILE = 'D:\Zichen_Projects\microtask_GEO\DO\HumanDO.obo'
	ontology_id_terms, name2id, G = obo_parser.obo2network(DO_FILE)
	return G.to_undirected()

def _plot_roc_similar_dz(id_pairs, G, spl_cutoff=1):
	## compute TPR, FPR, auroc for similar dz
	## G is the network for disease ontology
	## spl_cutoff shortest path length cutoff
	fp_list = []
	tp_list = []
	fp = 0
	tp = 0
	for uid_i, uid_j in id_pairs:
		gse_i, gse_j = d_uid_gse[uid_i], d_uid_gse[uid_j]
		if gse_i != gse_j: # exclude signatures from the same GEO study
			doid_i, doid_j = d_uid_doid[uid_i], d_uid_doid[uid_j]
			umls_i, umls_j = d_uid_umls[uid_i], d_uid_umls[uid_j]
			if doid_i == doid_j and None not in [doid_i, doid_j]:
				tp += 1
			elif umls_i == umls_j and None not in [umls_i, umls_j]:
				tp += 1
			elif G.has_node(doid_i) and G.has_node(doid_j):
				if nx.shortest_path_length(G, doid_i, doid_j) <= spl_cutoff:
					tp += 1
				else:
					fp +=1
			else:
				fp += 1
			tp_list.append(tp)
			fp_list.append(fp)
	max_TP = float(tp_list[-1])
	max_FP = float(fp_list[-1])
	TPRs = np.array(tp_list)/max_TP
	FPRs = np.array(fp_list)/max_FP
	auroc = auc(FPRs, TPRs)
	return auroc, FPRs, TPRs	

def _plot_roc_same_gene(id_pairs):
	## compute TPR, FPR, auroc for same gene
	fp_list = []
	tp_list = []
	fp = 0
	tp = 0
	for uid_i, uid_j in id_pairs:
		gse_i, gse_j = d_uid_gse[uid_i], d_uid_gse[uid_j]
		if gse_i != gse_j: # exclude signatures from the same GEO study
			doid_i, doid_j = d_uid_hs[uid_i], d_uid_hs[uid_j]
			umls_i, umls_j = d_uid_mm[uid_i], d_uid_mm[uid_j]
			if doid_i == doid_j and 'NULL' not in [doid_i, doid_j]:
				tp += 1
			elif umls_i == umls_j and 'NULL' not in [umls_i, umls_j]:
				tp += 1
			else:
				fp += 1
			tp_list.append(tp)
			fp_list.append(fp)
	max_TP = float(tp_list[-1])
	max_FP = float(fp_list[-1])
	TPRs = np.array(tp_list)/max_TP
	FPRs = np.array(fp_list)/max_FP
	auroc = auc(FPRs, TPRs)
	return auroc, FPRs, TPRs

def load_ppi():
	# 'D:\Zichen_Projects\\bias_correction\PPIN\PPIN_to_months_low.p',
	# 'D:\Zichen_Projects\\bias_correction\PPIN\PPIN_to_months.p',
	return pickle.load(open('D:\Zichen_Projects\\bias_correction\PPIN\PPIN_to_months_low.p', 'rb'))
	# return pickle.load(open('D:\Zichen_Projects\\bias_correction\PPIN\PPIN_to_months.p', 'rb'))

def _plot_roc_similar_gene(id_pairs, G, spl_cutoff=1):
	## compute TPR, FPR, auroc for similar dz
	## G is the PPI network
	## spl_cutoff shortest path length cutoff
	fp_list = []
	tp_list = []
	fp = 0
	tp = 0
	for uid_i, uid_j in id_pairs:
		gse_i, gse_j = d_uid_gse[uid_i], d_uid_gse[uid_j]
		if gse_i != gse_j: # exclude signatures from the same GEO study
			hs_i, hs_j = d_uid_hs[uid_i], d_uid_hs[uid_j]
			mm_i, mm_j = d_uid_mm[uid_i], d_uid_mm[uid_j]
			if hs_i == hs_j and None not in [hs_i, hs_j]:
				tp += 1
			elif mm_i == mm_j and None not in [mm_i, mm_j]:
				tp += 1
			elif G.has_node(hs_i) and G.has_node(hs_j):
				try:
					if nx.shortest_path_length(G, hs_i, hs_j) <= spl_cutoff:
						tp += 1
					else:
						fp +=1
				except nx.exception.NetworkXNoPath:
					fp += 1

			else:
				fp += 1
			tp_list.append(tp)
			fp_list.append(fp)
	max_TP = float(tp_list[-1])
	max_FP = float(fp_list[-1])
	TPRs = np.array(tp_list)/max_TP
	FPRs = np.array(fp_list)/max_FP
	auroc = auc(FPRs, TPRs)
	return auroc, FPRs, TPRs

def load_drug_similarities(cutoff):
	# fn = 'drug_smiles-morgan-dice.txt.gz'
	fn = 'drug_smiles-maccs-dice.txt.gz'
	# fn = 'drug_smiles-morgan_fb-tanimoto.txt.gz'
	# fn = 'drug_smiles-torsion-tanimoto.txt.gz'
	# fn = 'drug_smiles-rdkit-tanimoto.txt.gz'
	gs = nx.Graph()
	with gzip.open(fn) as f:
		for line in f:
			uid_i, uid_j, score = line.strip().split('\t')
			uid_i, uid_j = int(uid_i), int(uid_j)
			score = float(score)
			if score > cutoff:
				gs.add_edge(uid_i, uid_j)
	return gs

def _plot_roc_similar_drug(id_pairs, G):
	fp_list = []
	tp_list = []
	fp = 0
	tp = 0
	for uid_i, uid_j in id_pairs:
		gse_i, gse_j = d_uid_gse[uid_i], d_uid_gse[uid_j]
		if gse_i != gse_j: # exclude signatures from the same GEO study
			if G.has_edge(uid_i, uid_j):
				tp += 1
			else:
				fp += 1
		tp_list.append(tp)
		fp_list.append(fp)
	max_TP = float(tp_list[-1])
	max_FP = float(fp_list[-1])
	TPRs = np.array(tp_list)/max_TP
	FPRs = np.array(fp_list)/max_FP
	auroc = auc(FPRs, TPRs)
	return auroc, FPRs, TPRs

def plot_roc(signed_jaccard_fn, ax, absolute=False, plot=True, label=None, color=None):
	id_pairs = read_pairwise_from_file(signed_jaccard_fn, absolute=absolute)
	## diseases:
	# auroc, FPRs, TPRs = _plot_roc_same_dz(id_pairs)
	# G= load_do_network()
	# auroc, FPRs, TPRs = _plot_roc_similar_dz(id_pairs, G, spl_cutoff=3)
	## genes:
	# auroc, FPRs, TPRs = _plot_roc_same_gene(id_pairs)
	# G= load_ppi()
	# auroc, FPRs, TPRs = _plot_roc_similar_gene(id_pairs, G, spl_cutoff=2)
	## drugs:
	G = load_drug_similarities(0.9)
	auroc, FPRs, TPRs = _plot_roc_similar_drug(id_pairs, G)

	if plot:
		label += ", AUC = %.3f" % auroc
		print label
		ax.plot(FPRs, TPRs, label=label, color=color, lw=2)
		ax.set_xlabel('False Positive Rate',fontsize=20)
		ax.set_ylabel('True Positive Rate',fontsize=20)
		ax.legend(loc='lower right',prop={'size':16})
	return auroc

## plot roc curves for recovering known connections between signatures
'''
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
# plot_roc('signed_jaccard_839_dz_unique_entries.txt.gz', ax, absolute=False, label='signed jaccard', color='b')
# plot_roc('signed_jaccard_839_dz_unique_entries.txt.gz', ax, absolute=True, label='abs(signed jaccard)', color='r')

# plot_roc('signed_jaccard_2460_gene_unique_entries.txt.gz', ax, absolute=False, label='signed jaccard', color='b')
# plot_roc('signed_jaccard_2460_gene_unique_entries.txt.gz', ax, absolute=True, label='abs(signed jaccard)', color='r')

plot_roc('signed_jaccard_906_drug_unique_entries.txt.gz', ax, absolute=False, label='signed jaccard', color='b')
plot_roc('signed_jaccard_906_drug_unique_entries.txt.gz', ax, absolute=True, label='abs(signed jaccard)', color='r')


enlarge_tick_fontsize(ax, 16)
plt.show()
'''


## plot embedding for the adjacency matrix
'''
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
# sam.plot_embedding('TSNE')
# embedding = sam.tsne()
# embedding = np.loadtxt('tsne_signed_jaccard_4066_unique_entries.txt.gz.txt')
#### color by entry type
colors = []
categories = []
for prefix_id in unique_entries:
	if prefix_id.startswith('dz'):
		colors.append('r')
		categories.append('disease')
	elif prefix_id.startswith('drug'):
		colors.append('b')
		categories.append('drug')
	else:
		colors.append('k')
		categories.append('gene')

# fig = plt.figure(figsize=(10,10))
# ax = fig.add_subplot(111)
# ax.scatter(embedding[:,0], embedding[:,1], c=colors)
# plt.show()

clustergram(sam.to_csr_matrix().toarray(), row_groups=categories, col_groups=categories, display_range=0.05, 
	colorkey='signed jaccard',
	row_pdist='cosine', col_pdist='cosine',
	row_linkage='complete', col_linkage='complete'
	)

'''

