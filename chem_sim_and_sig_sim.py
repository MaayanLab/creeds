## use signature similarity to benchmark different chemical similarity measures
## to address issue#5
## created on 7/1/2015

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
from scipy.stats import pearsonr, spearmanr

from GEOentry import *
sys.path.append('C:\Users\Zichen\Documents\\bitbucket\maayanlab_utils')
from fileIO import file2list, mysqlTable2dict, write_df
from plots import enlarge_tick_fontsize, COLORS10

sys.path.append('C:\Users\Zichen\Documents\\bitbucket\\natural_products')
from SparseAdjacencyMat import SparseAdjacencyMat, sam_correlation


d_gds_gse = mysqlTable2dict('maaya0_crowdsourcing', 'geo2enrichr_dz', 0, 1)

## load cleaned drugs
global d_uid_gse
d_uid_geoid = mysqlTable2dict('maaya0_crowdsourcing', 'geo2enrichr_drug', -1, 0)

d_uid_gse = {} # convert gds to gse
for uid, geoid in d_uid_geoid.items():
	if geoid in d_gds_gse: gse = d_gds_gse[geoid]
	else: gse = geoid
	d_uid_gse[uid] = gse

global fp_types, metrics
fp_types = ['maccs', 'morgan', 'morgan_fb', 'rdkit', 'torsion']
metrics = ['tanimoto', 'dice']


def read_sam_from_file(fn, absolute=False, filter_same_gse=True, chem=True, cutoff=-2):
	## read signed jaccard file
	## return sam instance
	mat = {}
	with gzip.open(fn, 'rb') as f:
		for line in f:
			sl = line.strip().split('\t')
			if not chem: # $1 and $2 is prefix_id
				uid_i, uid_j = map(lambda x: int(x.split(':')[1]), sl[:2])
			else: # $1 and $2 is uid
				uid_i, uid_j = map(int, sl[:2])
			score = float(sl[2])
			if absolute:
				score = abs(score)

			if score > cutoff: # a cutoff for score, used when reading a gold standard
		
				if filter_same_gse: # to filter out signatures from the same GSE
					gse_i, gse_j = d_uid_gse[uid_i], d_uid_gse[uid_j]
					if gse_i != gse_j:
						mat[uid_i, uid_j] = score
				else:
					mat[uid_i, uid_j] = score
			
	fn = fn.split('/')[-1]
	return SparseAdjacencyMat(mat, fn)

def calc_all_corrs():
	# calculate all correlations
	sam_sig_abs = read_sam_from_file('signed_jaccard_906_drug_unique_entries.txt.gz', absolute=True, chem=False, filter_same_gse=False)
	sam_sig = read_sam_from_file('signed_jaccard_906_drug_unique_entries.txt.gz', absolute=False, chem=False, filter_same_gse=False)

	corr_mat = np.zeros((10, 4))
	i = 0

	header = ['chem_metrics', 'abs(signed jaccard) pearson', 'signed jaccard pearson', 
		'abs(signed jaccard) spearman', 'signed jaccard spearman']
	rids = []
	for fp_type in fp_types:
		for metric in metrics:
			chem_fn = 'drug_smiles-%s-%s.txt.gz' %(fp_type, metric)
			sam_chem = read_sam_from_file(chem_fn, filter_same_gse=False)
			corr_row = [sam_correlation(sam_chem, sam_sig_abs, type='pearson'),
						sam_correlation(sam_chem, sam_sig, type='pearson'),
						sam_correlation(sam_chem, sam_sig_abs, type='spearman'),
						sam_correlation(sam_chem, sam_sig, type='spearman'),
						]
			corr_row = np.array(corr_row)
			rids.append('-'.join([fp_type, metric]))
			print fp_type, metric, corr_row
			corr_mat[i] = corr_row
			i += 1
	write_df(corr_mat, rids, header, 'correlation_matrix_chemfps_drug_sigs.txt')
	return corr_mat


def plot_roc(sam, gs, ax, label=None, color=None, ls='-', plot=True):
	## plot a roc curve on ax using sam as input and a gs (also a sam) as gold standard
	## sort sam
	id_pairs = sam.d.keys()
	scores = sam.d.values()
	id_pairs = np.array(id_pairs)
	scores = np.array(scores)
	srt_idx = np.argsort(scores)[::-1] # large to small
	del sam, scores
	id_pairs = id_pairs[srt_idx]

	fp_list = []
	tp_list = []
	fp = 0
	tp = 0
	for uid_i, uid_j in id_pairs:
		if gs.has_edge(uid_i, uid_j):
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
	if plot:
		label += ", AUC = %.3f" % auroc
		print label
		ax.plot(FPRs, TPRs, label=label, color=color, lw=2)
		ax.set_xlabel('False Positive Rate',fontsize=20)
		ax.set_ylabel('True Positive Rate',fontsize=20)
		ax.legend(loc='lower right',prop={'size':16})
	return auroc

## calc all correlations
# calc_all_corrs()

## plot ROC curves
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)

gs = read_sam_from_file('signed_jaccard_906_drug_unique_entries.txt.gz', 
	absolute=False, chem=False, filter_same_gse=True, cutoff=0.08)
print len(gs)

i = 0
# lss = ['-', '--'] * 5
for fp_type in fp_types:
	# for metric in metrics:
	metric = 'tanimoto'
	chem_fn = 'drug_smiles-%s-%s.txt.gz' %(fp_type, metric)
	sam_chem = read_sam_from_file(chem_fn, filter_same_gse=True)
	plot_roc(sam_chem, gs, ax, label='%s' %(fp_type), color=COLORS10[i], ls='-', plot=True)
	i += 1

plt.show()
