## To plot a ROC curve for each dz/gene/drug entry 
## to show their ability to recover known targets(dz_genes/PPIs/targets)
## to address issue#1
## created on 6/25/2015

import os, sys
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 ## Output Type 3 (Type3) or Type 42 (TrueType)
rcParams['font.sans-serif'] = 'Arial'

import urllib2, json, MySQLdb
import cPickle as pickle
from itertools import combinations
from collections import Counter
from pprint import pprint
from sklearn.metrics import auc

from GEOentry import *
# from gene_convertion import *
from gene_convertion import _humanize
from jsons2gmt import get_valid_entries

sys.path.append('C:\Users\Zichen\Documents\\bitbucket\maayanlab_utils')
from fileIO import file2list, mysqlTable2dict, sqliteTable2dict, read_gmt
from plots import enlarge_tick_fontsize


def _plot_roc(sorted_list, gs):
	fp_list = []
	tp_list = []
	fp = 0
	tp = 0
	for gene in sorted_list:
		if gene in gs:
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

def _read_json_entry(fn,absolute=True):
	"""
	read a json file and return a sorted list of genes
	"""
	entry = json2entry(fn, meta_only=False)
	sorted_list = entry.to_sorted_gene_list(absolute=absolute)
	return map(_humanize, sorted_list) # humanize mouse genes

def make_gs(fn, val_idx, key_idx, header=False):
	# to store the gold standard from a file
	G = nx.Graph()
	with open (fn) as f:
		if header:
			next(f)
		for line in f:
			sl = line.strip().split('\t')
			val, key = sl[val_idx], sl[key_idx]
			G.add_edge(val, key)
	return G

def plot_roc(uid, key, G_gs, ax=None, label=None, color=None, plot=True):
	# key is the node name in G_gs
	json_fn = str(uid) + '.json'
	sorted_list = _read_json_entry(json_fn)
	gs = G_gs.neighbors(key)
	auroc, FPRs, TPRs = _plot_roc(sorted_list, gs)
	if plot:
		label += ", AUC = %.3f" % auroc
		print label
		ax.plot(FPRs, TPRs, label=label, color=color, lw=2)
		ax.set_xlabel('False Positive Rate',fontsize=20)
		ax.set_ylabel('True Positive Rate',fontsize=20)
		ax.legend(loc='lower right',prop={'size':16})
	return auroc

# valid_entries = get_valid_entries('output/microtask_gene_jsons', 'valid_gene_entries.pkl')
valid_entries = get_valid_entries('output/microtask_dz_jsons', 'valid_dz_entries.pkl')
# valid_entries = get_valid_entries('output/microtask_drug_jsons', 'valid_drug_entries.pkl')

# print valid_entries.keys()[:5]

## load dz id
d_uid_doid = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_dzs', 0, 1)
d_uid_umls = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_dzs', 0, 2)
# d_uid_name = mysqlTable2dict('maaya0_crowdsourcing', 'geo2enrichr_dz', -1, 3)
d_doid_name = dict(zip(file2list('D:\Zichen_Projects\microtask_GEO\DISEASES\human_disease_filtered_combined.tsv', 2),
	file2list('D:\Zichen_Projects\microtask_GEO\DISEASES\human_disease_filtered_combined.tsv', 3)))

def write_dz_genes_aurocs(gs_type, valid_entries, d_uid_doid):
	d_gs_type_fn = {
		'knowledge': 'D:\Zichen_Projects\microtask_GEO\DISEASES\human_disease_knowledge_filtered.tsv',
		'experiments': 'D:\Zichen_Projects\microtask_GEO\DISEASES\human_disease_experiments_filtered.tsv',
		'textmining': 'D:\Zichen_Projects\microtask_GEO\DISEASES\human_disease_textmining_filtered.tsv',
		'combined': 'D:\Zichen_Projects\microtask_GEO\DISEASES\human_disease_filtered_combined.tsv',
	}
	
	outfn = 'all_AUROCs_%s_dz_genes.txt' % gs_type
	G_k = make_gs(d_gs_type_fn[gs_type], 1, 2)
	with open (outfn, 'w') as out:
		for uid in valid_entries:
			doid = d_uid_doid[uid]
			if doid is not None and G_k.has_node(doid):
				if len(G_k.neighbors(doid)) > 10: # at least 10 known disease genes
					auroc = plot_roc(uid, doid, G_k, plot=False)
					print uid, doid, auroc, len(G_k.neighbors(doid))
					out.write('\t'.join(map(str, [uid, doid, auroc, len(G_k.neighbors(doid))])) + '\n')
	return

## write all AUROC for recovering dz_genes using signatures
# write_dz_genes_aurocs('knowledge', valid_entries, d_uid_doid)
# write_dz_genes_aurocs('combined', valid_entries, d_uid_doid)
# write_dz_genes_aurocs('experiments', valid_entries, d_uid_doid)
# write_dz_genes_aurocs('textmining', valid_entries, d_uid_doid)

## analyse, stratify AUROCs
import seaborn as sns
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)

### violinplot for all AUROCs
'''
all_aurocs = []

all_aurocs.append( map(float, file2list('all_AUROCs_experiments_dz_genes.txt', 2) ) )
all_aurocs.append( map(float, file2list('all_AUROCs_knowledge_dz_genes.txt', 2) ) )
all_aurocs.append( map(float, file2list('all_AUROCs_textmining_dz_genes.txt', 2) ) )
all_aurocs.append( map(float, file2list('all_AUROCs_combined_dz_genes.txt', 2) ) )

sns.violinplot(all_aurocs, names=['experiments', 'knowledge', 'textmining', 'combined'], ax=ax)
enlarge_tick_fontsize(ax,16)
ax.set_ylabel('AUROC', fontsize=18)
plt.show()
'''

### scatter plot for AUROC vs number of known dz_genes
# ax.scatter(map(float, file2list('all_AUROCs_combined_dz_genes.txt', 2) ), 
# 	map(float, file2list('all_AUROCs_combined_dz_genes.txt', 3) ) )
# ax.set_xlabel('AUROCs', fontsize=18)
# ax.set_ylabel('Number of known disease genes', fontsize=18)
# enlarge_tick_fontsize(ax,16)
# plt.show()

### max for individual diseases
def read_auroc_by_dz(fn, d_uid_doid):
	# read the file output by 'write_dz_genes_aurocs'
	# and store the aurocs for each dz
	d_doid_aurocs = {}
	with open(fn) as f:
		for line in f:
			sl = line.strip().split('\t')
			uid = int(sl[0])
			auroc = float(sl[2])
			doid = d_uid_doid[uid]
			if doid not in d_doid_aurocs:
				d_doid_aurocs[doid] = []
			d_doid_aurocs[doid].append(auroc)
	return d_doid_aurocs

d_doid_aurocs = read_auroc_by_dz('all_AUROCs_combined_dz_genes.txt', d_uid_doid)
max_aurocs_per_dzs = np.array([max(aurocs) for aurocs in d_doid_aurocs.values()])
doids = np.array(d_doid_aurocs.keys())
dz_names = np.array([d_doid_name[doid] for doid in doids])

d_doid_num_dz_genes = dict(zip( file2list('all_AUROCs_combined_dz_genes.txt', 1), 
	map(int, file2list('all_AUROCs_combined_dz_genes.txt', 3) ) ))

num_dz_genes = np.array([d_doid_num_dz_genes[doid] for doid in doids])

ax.scatter(max_aurocs_per_dzs, num_dz_genes)
ax.set_xlabel('max AUROCs', fontsize=18)
ax.set_ylabel('Number of known disease genes', fontsize=18)
enlarge_tick_fontsize(ax,16)
plt.show()


### plot histogram for max AUROCs for each dz
# sns.distplot(max_aurocs_per_dzs, rug=True, hist=False, ax=ax)
# ax.set_xlabel('max AUROC for each disease', fontsize=18)
# enlarge_tick_fontsize(ax,16)
# plt.show()


# print max_aurocs_per_dzs[max_aurocs_per_dzs > 0.7]
print doids[max_aurocs_per_dzs>0.7]
print dz_names[max_aurocs_per_dzs>0.7]

# print max_aurocs_per_dzs[max_aurocs_per_dzs > 0.5]
print doids[max_aurocs_per_dzs<0.45]
print dz_names[max_aurocs_per_dzs<0.45]
