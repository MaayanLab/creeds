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
from pymongo import MongoClient

from GEOentry import *
sys.path.append('C:\Users\Zichen\Documents\\bitbucket\maayanlab_utils')
from fileIO import file2list, mysqlTable2dict, write_df
from plots import enlarge_tick_fontsize, COLORS10
import ROC_test

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
# metrics = ['tanimoto', 'dice']
metrics = ['tanimoto']

def get_uid_category(prefix, field):
	# get d_uid_category from mongodb
	d_uid_category = {}
	client = MongoClient('mongodb://127.0.0.1:27017/')
	# client = MongoClient('mongodb://146.203.54.131:27017/')
	db = client['microtask_signatures']
	COLL = db['signatures']
	all_uids = [uid for uid in COLL.distinct('id') if uid.startswith(prefix)]
	projection = {'_id':False, 'chdir':False,'limma':False, 'fold_changes':False}
	for uid in all_uids:
		doc = COLL.find_one({'id': uid}, projection)
		id = int(uid.split(':')[1])
		d_uid_category[id] = doc[field]
	return d_uid_category


def read_sam_from_file(fn, absolute=False, filter_same_gse=True, chem=True, cutoff=-2, d_uid_category=None, keep_same=True):
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
						if d_uid_category is not None: # keep only the save category
							try:
								cat_i, cat_j = d_uid_category[uid_i], d_uid_category[uid_j]
							except:
								cat_i, cat_j = 0, 1
							if keep_same and cat_i == cat_j: # only keep uid pairs of the save category
								mat[uid_i, uid_j] = score
							elif (not keep_same) and cat_i != cat_j:
								mat[uid_i, uid_j] = score
						else:
							mat[uid_i, uid_j] = score
				else:
					mat[uid_i, uid_j] = score
			
	fn = fn.split('/')[-1]
	return SparseAdjacencyMat(mat, fn)

def calc_all_corrs():
	# calculate all correlations
	d_uid_platform = get_uid_category('drug', 'platform')

	# sam_sig_abs = read_sam_from_file('signed_jaccard_906_drug_unique_entries.txt.gz', absolute=True, chem=False
	# 	, d_uid_category=d_uid_platform, keep_same=False)
	sam_sig_diff = read_sam_from_file('signed_jaccard_906_drug_unique_entries.txt.gz', absolute=False, chem=False,
		filter_same_gse=True, d_uid_category=d_uid_platform, keep_same=False)
	sam_sig_same = read_sam_from_file('signed_jaccard_906_drug_unique_entries.txt.gz', absolute=False, chem=False,
		filter_same_gse=True, d_uid_category=d_uid_platform, keep_same=True)

	corr_mat = np.zeros((10, 4))
	i = 0

	header = ['chem_metrics', 'abs(signed jaccard) pearson', 'signed jaccard pearson', 
		'abs(signed jaccard) spearman', 'signed jaccard spearman']
	rids = []
	for fp_type in fp_types:
		for metric in metrics:
			chem_fn = 'drug_smiles-%s-%s.txt.gz' %(fp_type, metric)
			sam_chem = read_sam_from_file(chem_fn, filter_same_gse=True)
			# corr_row = [sam_correlation(sam_chem, sam_sig_abs, type='pearson'),
			# 			sam_correlation(sam_chem, sam_sig, type='pearson'),
			# 			sam_correlation(sam_chem, sam_sig_abs, type='spearman'),
			# 			sam_correlation(sam_chem, sam_sig, type='spearman'),
			# 			]
			# corr_row = np.array(corr_row)
			# rids.append('-'.join([fp_type, metric]))
			print fp_type, metric#, corr_row
			corr, vals1, vals2 = sam_correlation(sam_chem, sam_sig_diff, type='pearson', return_arrays=True)
			plt.scatter(vals1, vals2, color='b', alpha=0.1, label='diff platform, corr=%.3f'%corr)
			print corr
			corr, vals1, vals2 = sam_correlation(sam_chem, sam_sig_same, type='pearson', return_arrays=True)
			plt.scatter(vals1, vals2, color='r', alpha=0.1, label='same platform, corr=%.3f'%corr)
			plt.ylim([-0.4,1])
			print corr
			plt.legend()
			plt.show()
	# 		corr_mat[i] = corr_row
	# 		i += 1
	# write_df(corr_mat, rids, header, 'correlation_matrix_chemfps_drug_sigs_same_organism.txt')
	return corr_mat

def ge_chem_corr_plot(category, fp_type, ax, sam_sig_diff, sam_sig_same):
	## scatter plot to show the correlation between GE and chem scores 
	## divided by in the same/diff category
	metric = 'tanimoto'

	chem_fn = 'drug_smiles-%s-%s.txt.gz' %(fp_type, metric)
	sam_chem = read_sam_from_file(chem_fn, filter_same_gse=True)
	corr, vals1, vals2 = sam_correlation(sam_chem, sam_sig_diff, type='pearson', return_arrays=True)
	ax.scatter(vals1, vals2, color='b', alpha=0.2, label='corr=%.3f'%(corr), rasterized=True)
	corr, vals1, vals2 = sam_correlation(sam_chem, sam_sig_same, type='pearson', return_arrays=True)
	ax.scatter(vals1, vals2, color='r', alpha=0.2, label='corr=%.3f'%(corr), rasterized=True)
	ax.legend(prop={"size":16},loc='upper left')
	ax.set_ylim([-0.4,1])
	ax.set_xlim([-0.1,1.1])
	return

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

def plot_roc2(sam, gs):
	## plot a roc curve on ax using sam as input and a gs (also a sam) as gold standard
	## sort sam
	## and return a list of response [0,1] and predictor for roc_test
	id_pairs = sam.d.keys()
	scores = sam.d.values()

	response = []
	for uid_i, uid_j in id_pairs:
		if gs.has_edge(uid_i, uid_j):
			response.append(1)
		else:
			response.append(0)

	return response, scores


## calc all correlations
# calc_all_corrs()
## scatter plot for same/diff categories
# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
'''
fig, ax_mat = plt.subplots(2, 5, sharex='col', sharey='row', figsize=(20,8))

fp_types = {'maccs':'MACCS', 'morgan': 'ECFP4', 'morgan_fb': 'ECFP4-FB', 'rdkit': 'RDKit', 'torsion':'Torsions'}

category = 'platform'
d_uid_platform = get_uid_category('drug', category)
sam_sig_diff = read_sam_from_file('signed_jaccard_906_drug_unique_entries.txt.gz', absolute=False, chem=False,
	filter_same_gse=True, d_uid_category=d_uid_platform, keep_same=False)
sam_sig_same = read_sam_from_file('signed_jaccard_906_drug_unique_entries.txt.gz', absolute=False, chem=False,
	filter_same_gse=True, d_uid_category=d_uid_platform, keep_same=True)

for i, fp_type in enumerate(['maccs', 'morgan', 'morgan_fb', 'rdkit', 'torsion']):
	ax = ax_mat[0, i]
	ge_chem_corr_plot(category, fp_type, ax, sam_sig_diff, sam_sig_same)
	fp_show = fp_types[fp_type]
	ax.set_title(fp_show, fontsize=20)

category = 'organism'
d_uid_platform = get_uid_category('drug', category)
sam_sig_diff = read_sam_from_file('signed_jaccard_906_drug_unique_entries.txt.gz', absolute=False, chem=False,
	filter_same_gse=True, d_uid_category=d_uid_platform, keep_same=False)
sam_sig_same = read_sam_from_file('signed_jaccard_906_drug_unique_entries.txt.gz', absolute=False, chem=False,
	filter_same_gse=True, d_uid_category=d_uid_platform, keep_same=True)

for i, fp_type in enumerate(['maccs', 'morgan', 'morgan_fb', 'rdkit', 'torsion']):
	ax = ax_mat[1, i]
	ge_chem_corr_plot(category, fp_type, ax, sam_sig_diff, sam_sig_same)
	ax.set_xlabel('Tanimoto coefficient', fontsize=18)

ax_mat[0,0].set_ylabel('Signed Jaccard index', fontsize=18)
ax_mat[1,0].set_ylabel('Signed Jaccard index', fontsize=18)

ax_mat[0,0].set_title(fp_show, fontsize=20)

fig.tight_layout()
plt.show()
'''


# '''
## plot ROC curves
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)

d_uid_platform = get_uid_category('drug', 'platform')

# gs = read_sam_from_file('signed_jaccard_906_drug_unique_entries.txt.gz', 
# 	absolute=False, chem=False, filter_same_gse=True, cutoff=0.1)

gs = read_sam_from_file('signed_jaccard_906_drug_unique_entries.txt.gz', 
	absolute=False, chem=False, filter_same_gse=True, cutoff=0.1,
	d_uid_category=d_uid_platform, keep_same=False)


print len(gs)

i = 0
# lss = ['-', '--'] * 5
all_predictions = []
for fp_type in fp_types:
	# for metric in metrics:
	metric = 'tanimoto'
	chem_fn = 'drug_smiles-%s-%s.txt.gz' %(fp_type, metric)
	sam_chem = read_sam_from_file(chem_fn, filter_same_gse=True)
	# plot_roc(sam_chem, gs, ax, label='%s' %(fp_type), color=COLORS10[i], ls='-', plot=True)
	response, prediction = plot_roc2(sam_chem, gs)
	all_predictions.append(prediction)
	i += 1


## perform DeLong's test on AUC
for i,j in combinations(range(len(fp_types)), 2):
	pval, auc1, auc2 = ROC_test.roc_test(response, all_predictions[i], all_predictions[j])
	print fp_types[i], auc1, fp_types[j], auc2, pval


plt.show()
# '''
