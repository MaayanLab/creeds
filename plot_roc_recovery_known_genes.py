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
from plots import enlarge_tick_fontsize, COLORS10
from plot_DRR import plt_DRR, get_DRR_plt_auc


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
	# return map(_humanize, sorted_list) # humanize mouse genes
	return map(lambda x:x.upper(), sorted_list)

def _read_limma_from_json(uid):
	prefix, id = uid.split(':')
	fn = 'C://Users/Zichen/Documents/bitbucket/microtask_GEO/output/microtask_%s_jsons_limma_norm/%s.json' %(prefix, id)
	data = json.load(open(fn, 'rb'))
	return map(lambda x : _humanize(x[0]), data['limma'])
	

from pymongo import MongoClient
def _read_limma(uid):
	client = MongoClient('mongodb://146.203.54.131:27017/')
	db = client['microtask_signatures']
	COLL = db['signatures']
	projection = {'_id':False, 'limma':True}
	doc = COLL.find_one({'id':uid}, projection)
	if doc is not None:
		if 'limma' in doc:
			return map(lambda x: x.upper(), doc['limma']['genes'])


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
	# sorted_list = _read_json_entry(json_fn)
	sorted_list = _read_limma('dz:' + str(uid))
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

def _get_scaled_ranks(sorted_list, gs):
	# alternative to _plot_roc function to benchmark a sorted gene list
	ranks = np.arange(1, len(sorted_list) + 1)
	ranks = ranks / float(ranks[-1]) # scaled ranks
	mask = np.in1d(sorted_list, list(gs))
	return ranks[mask].tolist()	

def get_scaled_ranks(uid, key, G_gs, signature_type="limma"):
	## wrapper
	# key is the node name in G_gs
	# sorted_list = _read_limma('dz:' + str(uid))
	if signature_type == 'limma':
		# sorted_list = _read_limma_from_json('dz:' + str(uid))
		sorted_list = _read_limma('dz:' + str(uid))
	elif signature_type == 'chdir':
		json_fn = str(uid) + '.json'
		sorted_list = _read_json_entry(json_fn)
	if type(G_gs) == dict:
		gs = G_gs[key]
	else:
		gs = G_gs.neighbors(key)
	sr = _get_scaled_ranks(sorted_list, gs)
	return sr


# valid_entries = get_valid_entries('output/microtask_gene_jsons', 'valid_gene_entries.pkl')
valid_entries = get_valid_entries('output/microtask_dz_jsons', 'valid_dz_entries.pkl')
# valid_entries = get_valid_entries('output/microtask_drug_jsons', 'valid_drug_entries.pkl')

# print valid_entries.keys()[:5]

## load dz id
d_uid_doid = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_dzs', 0, 2)
d_uid_umls = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_dzs', 0, 3)
# d_uid_name = mysqlTable2dict('maaya0_crowdsourcing', 'geo2enrichr_dz', -1, 3)
d_doid_name = dict(zip(file2list('D:\Zichen_Projects\microtask_GEO\DISEASES\human_disease_filtered_combined.tsv', 2),
	file2list('D:\Zichen_Projects\microtask_GEO\DISEASES\human_disease_filtered_combined.tsv', 3)))

## load tf 
d_uid_TF = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_genes', 0, 1)

def write_dz_genes_aurocs(gs_type, valid_entries, d_uid_doid):
	d_gs_type_fn = {
		'knowledge': 'D:\Zichen_Projects\microtask_GEO\DISEASES\human_disease_knowledge_filtered.tsv',
		'experiments': 'D:\Zichen_Projects\microtask_GEO\DISEASES\human_disease_experiments_filtered.tsv',
		'textmining': 'D:\Zichen_Projects\microtask_GEO\DISEASES\human_disease_textmining_filtered.tsv',
		'combined': 'D:\Zichen_Projects\microtask_GEO\DISEASES\human_disease_filtered_combined.tsv',
	}
	
	outfn = 'all_AUROCs_limma_%s_dz_genes.txt' % gs_type
	G_k = make_gs(d_gs_type_fn[gs_type], 1, 2)
	with open (outfn, 'w') as out:
		for uid in valid_entries:
			doid = d_uid_doid[uid]
			if doid is not None and G_k.has_node(doid):
				if len(G_k.neighbors(doid)) > 10: # at least 10 known disease genes
					try:
						auroc = plot_roc(uid, doid, G_k, plot=False)
						print uid, doid, auroc, len(G_k.neighbors(doid))
						out.write('\t'.join(map(str, [uid, doid, auroc, len(G_k.neighbors(doid))])) + '\n')
					except:
						pass
	return

def plot_dz_genes_DRR(gs_type, valid_entries, d_uid_doid, ax, ax_right, signature_type=None, color=None):
	d_gs_type_fn = {
		'knowledge': 'D:\Zichen_Projects\microtask_GEO\DISEASES\human_disease_knowledge_filtered.tsv',
		'experiments': 'D:\Zichen_Projects\microtask_GEO\DISEASES\human_disease_experiments_filtered.tsv',
		'textmining': 'D:\Zichen_Projects\microtask_GEO\DISEASES\human_disease_textmining_filtered.tsv',
		'combined': 'D:\Zichen_Projects\microtask_GEO\DISEASES\human_disease_filtered_combined.tsv',
	}
	
	G_k = make_gs(d_gs_type_fn[gs_type], 1, 2)

	srs = []
	for uid in valid_entries:
		doid = d_uid_doid[uid]
		if doid is not None and G_k.has_node(doid):
			if len(G_k.neighbors(doid)) > 0: # at least 10 known disease genes
				print uid
				try:
					sr = get_scaled_ranks(uid, doid, G_k, signature_type=signature_type)
					srs.extend(sr)
				except:
					pass
	plt_DRR(srs, ax, ax_right, color=color, ls='-', label=signature_type)
	return

def plot_tf_genes_DRR(gmt, valid_entries, d_uid_TF, ax, ax_right, signature_type=None, color=None):
	## try to repeat limma beat chdir on kevin's TF set
	HOME = 'C:/Users/Zichen/'
	if gmt == 'ChEA':
		d_gmt = read_gmt(HOME+'/Documents/GitLab/gmt_enrichr/ChEA.gmt')
		sep = '-'
	elif gmt == 'ENCODE':
		d_gmt = read_gmt(HOME+'/Documents/GitLab/gmt_enrichr/ENCODE_TF_ChIP-seq_2015.gmt')
		sep = '_'
	elif gmt == 'TF-PPI':
		d_gmt = read_gmt(HOME+'/Documents/GitLab/gmt_enrichr/Transcription_Factor_PPIs.gmt')
		sep = ''
	d_tf_targets = {}
	for term, genes in d_gmt.items():
		tf = term.split(sep)[0] # ChEA
		if tf not in d_tf_targets:
			d_tf_targets[tf] = set(genes)
		else:
			d_tf_targets[tf] = set(genes) & d_tf_targets[tf]

	srs = []
	for uid in valid_entries:
		tf = d_uid_TF[uid]
		if tf is not None and tf in d_tf_targets:
			print uid, tf
			try:
				sr = get_scaled_ranks(uid, tf, d_tf_targets, signature_type=signature_type)
				srs.extend(sr)
				print len(srs)
			except:
				pass

	plt_DRR(srs, ax, ax_right, color=color, ls='-', label=signature_type)
	return

## write all AUROC for recovering dz_genes using signatures
# write_dz_genes_aurocs('knowledge', valid_entries, d_uid_doid)
# write_dz_genes_aurocs('combined', valid_entries, d_uid_doid)
# write_dz_genes_aurocs('experiments', valid_entries, d_uid_doid)
# write_dz_genes_aurocs('textmining', valid_entries, d_uid_doid)

## plot DRR 
fig = plt.figure(figsize=(9,9))
ax = fig.add_subplot(111)
ax_right = ax.twinx()

gs_type = 'knowledge'
# plot_dz_genes_DRR(gs_type, valid_entries, d_uid_doid, ax, ax_right, signature_type='chdir', color=COLORS10[0])
plot_dz_genes_DRR(gs_type, valid_entries, d_uid_doid, ax, ax_right, signature_type='limma', color=COLORS10[1])
# # plot_dz_genes_DRR('combined', valid_entries, d_uid_doid)

# ax.set_title(gs_type, fontsize=18)

# gmt = 'ChEA'
# d_uid_curator = mysqlTable2dict('maaya0_crowdsourcing', 'geo2enrichr', -1, -3)
# valid_entries = []

# for uid, curator in d_uid_curator.items():
# 	# if curator == 'Kevin' :
# 	valid_entries.append(uid)
# print len(valid_entries)

# plot_tf_genes_DRR(gmt,valid_entries, d_uid_TF, ax, ax_right, signature_type='chdir', color=COLORS10[0])
# plot_tf_genes_DRR(gmt,valid_entries, d_uid_TF, ax, ax_right, signature_type='limma', color=COLORS10[1])

# ax.set_title(gmt, fontsize=18)
ax.legend(loc='best')
plt.show()

## analyse, stratify AUROCs
'''
import seaborn as sns
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)

### violinplot for all AUROCs

# all_aurocs = []

# all_aurocs.append( map(float, file2list('all_AUROCs_limma_experiments_dz_genes.txt', 2) ) )
# all_aurocs.append( map(float, file2list('all_AUROCs_limma_knowledge_dz_genes.txt', 2) ) )
# all_aurocs.append( map(float, file2list('all_AUROCs_limma_textmining_dz_genes.txt', 2) ) )
# all_aurocs.append( map(float, file2list('all_AUROCs_limma_combined_dz_genes.txt', 2) ) )

# sns.violinplot(all_aurocs, names=['experiments', 'knowledge', 'textmining', 'combined'], ax=ax)
# enlarge_tick_fontsize(ax,16)
# ax.set_ylabel('AUROC', fontsize=18)
# plt.show()

### scatter plot for AUROC vs number of known dz_genes
# ax.scatter(map(float, file2list('all_AUROCs_limma_combined_dz_genes.txt', 2) ), 
# 	map(float, file2list('all_AUROCs_limma_combined_dz_genes.txt', 3) ) )
# ax.set_xlabel('AUROCs', fontsize=18)
# ax.set_ylabel('Number of known disease genes', fontsize=18)
# enlarge_tick_fontsize(ax,16)
# plt.show()

### scatter plot comparing limma and chdir

d_uid_auc_chdir = dict(zip(file2list('all_AUROCs_knowledge_dz_genes.txt',0), file2list('all_AUROCs_knowledge_dz_genes.txt',2)))
d_uid_auc_limma = dict(zip(file2list('all_AUROCs_limma_knowledge_dz_genes.txt',0), file2list('all_AUROCs_limma_knowledge_dz_genes.txt',2)))
d_uid_num_dz_genes = dict(zip(file2list('all_AUROCs_limma_knowledge_dz_genes.txt',0), file2list('all_AUROCs_limma_knowledge_dz_genes.txt',3))) 
shared_uids = set(d_uid_auc_chdir.keys()) & set(d_uid_auc_limma.keys())
ax.scatter([float(d_uid_auc_chdir[uid]) for uid in shared_uids], [float(d_uid_auc_limma[uid]) for uid in shared_uids],
	s=[int(d_uid_num_dz_genes[uid]) for uid in shared_uids], alpha=0.6)

ax.plot([0.2,1],[0.2,1], ls='--', color='grey')
ax.set_xlim([0.2,1])
ax.set_ylim([0.2,1])
ax.set_xlabel('AUROCs chdir', fontsize=18)
ax.set_ylabel('AUROCs limma', fontsize=18)
enlarge_tick_fontsize(ax,16)
plt.show()

'''

'''

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
'''