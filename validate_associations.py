import os, sys
import numpy as np
import cPickle as pickle
import networkx as nx
import matplotlib.pyplot as plt
from sklearn.metrics import auc
from scipy.stats import rankdata
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 ## Output Type 3 (Type3) or Type 42 (TrueType)
rcParams['font.sans-serif'] = 'Arial'

import seaborn as sns

import json
from GEOentry import *
sys.path.append('C:\Users\Zichen\Documents\\bitbucket\maayanlab_utils')
from fileIO import file2list
from plots import enlarge_tick_fontsize, COLORS10
from homologene import *


## parse HumanDO.obo to get the hierarchy of DOIDs
import obo_parser
DO_FILE = 'D:\Zichen_Projects\microtask_GEO\DO\HumanDO.obo'
_, _, G_DO = obo_parser.obo2network(DO_FILE)

## Gold standards: do_rif, DisGeNET, DISEASES
d_geneID_geneSymbol = id2symbol('9606')

def make_gs(fn, gene_idx, dz_idx, header=False):
	# to store the gold standard from a file
	G = nx.Graph()
	with open (fn) as f:
		if header:
			next(f)
		for line in f:
			sl = line.strip().split('\t')
			gene, dz = sl[gene_idx], sl[dz_idx]
			# gene, dz = sl[gene_idx], sl[dz_idx].split(':')[1]
			G.add_edge(gene, dz)
	return G

def prepare_predictions(edgelist_fn, d1, d2, G_DO=G_DO):
	# backtrack DOIDs to all parent DOIDs
	root = nx.topological_sort(G_DO)[0]
	all_edges = []
	with open (edgelist_fn) as f:
		for line in f:
			sl = line.strip().split('\t')
			gene = d1[int(sl[0])]
			dz = d2[int(sl[1])]

			if gene != None and dz != None:
				score = sl[2]
				# if G_DO.has_node(dz):
				# 	path = nx.shortest_path(G_DO, source=root, target=dz)
				# 	edges = [(gene, d, score) for d in path]
				# 	all_edges.extend(edges)
				# else: # not DOID
				all_edges.append((gene, dz, score))

	all_edges = np.array(all_edges)
	return all_edges[:, 0:2], all_edges[:, 2].astype(float)

def _validation(edges, G_gs):
	# validate overall gene-dz associations
	fps = []
	tps = []
	fp = 0
	tp = 0
	i = 0
	mask = []
	for gene, dz in edges:
		if G_gs.has_node(gene) and G_gs.has_node(dz):
			if G_gs.has_edge(gene,dz):
				tp += 1
			else:
				fp += 1
			tps.append(tp)
			fps.append(fp)
			mask.append(i)
		else:
			continue
		i += 1	
	max_tp = float(tps[-1])
	max_fp = float(fps[-1])
	# print 'max_tp=', max_tp
	return np.array(fps), np.array(tps), mask

def _validataion_per_dz(edges, G_gs):
	# validate gene-dz associations on dz-basis
	all_dzs = np.unique(edges[:, 1])
	d_dz_auc = {}
	for dz in set(all_dzs) & set(G_gs.nodes()):
		mask = np.in1d(edges[:, 1], [dz])
		fps, tps, _ = _validation(edges[mask], G_gs)
		try:
			auc_val = auc(fps/float(fps[-1]), tps/float(tps[-1]))
			d_dz_auc[dz] = auc_val
		except:
			pass
	return d_dz_auc



def _validation2(edges, ppi_fn):
	G = pickle.load(open(ppi_fn, 'rb'))
	fps = []
	tps = []
	fp = 0
	tp = 0
	for g1, g2 in edges:
		if g1 == g2:
			tp += 1
		else:
			if G.has_edge(g1,g2):
				tp += 1
			else:
				try:
					path_length = nx.shortest_path_length(G, source=g1, target=g2)
					if path_length == 2:
						tp += 1
					else:
						fp += 1
				except:
					fp += 1

		tps.append(tp)
		fps.append(fp)
	max_tp = float(tps[-1])
	max_fp = float(fps[-1])
	print 'max_tp=', max_tp
	return np.array(fps), np.array(tps)
	

def plot_roc(edgelist_fn, d1, d2, G_gs, ax, ls='-', label=None):
	## d1: d_uid_name for col1
	## d2: d_uid_name for col2

	edges, scores = prepare_predictions(edgelist_fn, d1, d2)

	fps1, tps1, mask = _validation(edges, G_gs)

	scores = scores[mask]
	edges = edges[mask]

	srt_idx = abs(scores).argsort()[::-1]
	# ranks1 = rankdata(scores, method='min')[::-1]
	# FPR1, TPR1 = ranks1/max(ranks1), tps1/float(tps1[-1])
	FPR1, TPR1 = fps1/float(fps1[-1]), tps1/float(tps1[-1])
	roc_auc1 = auc(FPR1, TPR1)
	# roc_auc1 = 0
	ax.plot(FPR1, TPR1, color=COLORS10[0], label='signed Jaccard' +' (AUC = %0.3f)' % roc_auc1, linewidth=2, ls=ls)


	edges, scores = prepare_predictions(edgelist_fn, d1, d2)
	srt_idx = abs(scores).argsort()[::-1]
	fps2, tps2, mask = _validation(edges[srt_idx], G_gs)

	scores = scores[mask]
	edges = edges[mask]

	srt_idx = abs(scores).argsort()[::-1]

	abs_scores_sorted = abs(scores[srt_idx])
	# ranks2 = rankdata(abs_scores_sorted, method='min')[::-1]
	# FPR2, TPR2 = ranks2/max(ranks2), tps2/float(tps2[-1])
	FPR2, TPR2 = fps2/float(fps2[-1]), tps2/float(tps2[-1])
	roc_auc2 = auc(FPR2, TPR2)
	# roc_auc2 = 0
	ax.plot(FPR2, TPR2, color=COLORS10[1], label='abs(signed Jaccard)'+' (AUC = %0.3f)' % roc_auc2, linewidth=2, ls=ls)
	
	ax.set_xlabel('False Positive Rate',fontsize=16)
	# ax.set_xlabel('signed Jaccard',fontsize=16)
	ax.set_ylabel('True Positive Rate',fontsize=16)
	# ax.set_xlabel('ranks',fontsize=16)
	# ax.set_ylabel('# True gene-dz associations',fontsize=16)
	# ax.set_xlim([-0.1, 0.1])
	# ax.set_ylim([-0.1, 0.1])
	enlarge_tick_fontsize(ax, 16)
	ax.legend(loc='lower right',prop={'size':14})

def box_auc(edgelist_fn, d1, d2, G_gs, ax):
	# box plot AUCs for each dz
	edges, scores = prepare_predictions(edgelist_fn, d1, d2)
	d_dz_auc1 = _validataion_per_dz(edges, G_gs)

	srt_idx = abs(scores).argsort()[::-1]
	d_dz_auc2 = _validataion_per_dz(edges[srt_idx], G_gs)

	np.random.shuffle(edges)
	d_dz_auc_s = _validataion_per_dz(edges, G_gs)

	sns.violinplot([d_dz_auc1.values(), d_dz_auc2.values(), d_dz_auc_s.values()], names=['signed Jaccard', 'abs(signed Jaccard)', 'shuffled'], ax=ax)
	ax.set_ylabel('AUROC per disease', fontsize=16)
	enlarge_tick_fontsize(ax, 14)
	print np.mean(d_dz_auc1.values()), np.mean(d_dz_auc2.values()), np.mean(d_dz_auc_s.values())
	return

## load dz id
d_uid_doid = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_dzs', 0, 1)
d_uid_umls = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_dzs', 0, 2)
# d_umls_uid = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_dzs', 2, 1)
# d_uid_dzname = mysqlTable2dict('maaya0_crowdsourcing', 'geo2enrichr_dz', -1, 3)

d_uid_gene = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_genes', 0, 1)

os.chdir('D:\Zichen_Projects\microtask_GEO\DISEASES')
G_k = make_gs('human_disease_knowledge_filtered.tsv', 1, 2)
G_k_full = make_gs('human_disease_knowledge_full.tsv', 1, 2)
# G_k = make_gs('human_disease_experiments_filtered.tsv', 1, 2)
# G_k_full = make_gs('human_disease_experiments_full.tsv', 1, 2)
# G_k = make_gs('human_disease_textmining_filtered.tsv', 1, 2)
# G_k_full = make_gs('human_disease_textmining_full.tsv', 1, 2)
# G_k = make_gs('human_disease_filtered_combined.tsv', 1, 2)
# G_k_full = make_gs('human_disease_full_combined.tsv', 1, 2)


# os.chdir('D:\Zichen_Projects\microtask_GEO\DisGeNET')
# G_k = make_gs('curated_gene_disease_associations.txt', 1, 3, header=True)
# G_k_full = make_gs('befree_gene_disease_associations.txt', 1, 3, header=True)


os.chdir('C:\\Users\Zichen\Documents\\bitbucket\\microtask_GEO\\output\\analyses_output')

PPIN_fns = [
'D:\Zichen_Projects\\bias_correction\PPIN\PPIN_to_months_low.p',
'D:\Zichen_Projects\\bias_correction\PPIN\PPIN_to_months.p',
]

fig = plt.figure(figsize=(12,6))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

# plot_roc('gene_umls_edgelist_signed_jaccard.txt', d_umls_gene , ax)

# plot_roc('gene_gene_edgelist_signed_jaccard.txt', PPIN_fns[0] , ax, ls='-' )
# plot_roc('gene_gene_edgelist_signed_jaccard.txt', PPIN_fns[1] , ax, ls='--' )

# plot_roc('gene_dz_edgelist_signed_jaccard.txt', d_uid_gene, d_uid_umls, G_k, ax1)
# plot_roc('gene_dz_edgelist_signed_jaccard.txt', d_uid_gene, d_uid_umls, G_k_full, ax2)
# box_auc('gene_dz_edgelist_signed_jaccard.txt', d_uid_gene, d_uid_doid, G_k, ax1)
# box_auc('gene_dz_edgelist_signed_jaccard.txt', d_uid_gene, d_uid_doid, G_k_full, ax2)


# dz_uids = map(int, set(file2list('gene_dz_edgelist_signed_jaccard.txt', 1)))
# dz_umls = [d_uid_umls[uid] for uid in dz_uids if uid in d_uid_umls]
# from collections import Counter
# d_umls_freq = dict(Counter(dz_umls))
# # from pprint import pprint
# # pprint(d_umls_freq)
# print sorted(d_umls_freq.items(), key=lambda x:x[1], reverse=True)[0:10]
# print len(d_umls_freq)



ax1.set_title('DISEASES knowledge filtered', fontsize=20)
ax2.set_title('DISEASES knowledge full', fontsize=20)
# ax1.set_title('DisGeNET curated', fontsize=20)
# ax2.set_title('DisGeNET befree', fontsize=20)


# fig.tight_layout()
# plt.show()

edges, scores = prepare_predictions('gene_dz_edgelist_signed_jaccard.txt', d_uid_gene, d_uid_doid)
srt_idx = abs(scores).argsort()[::-1]
d_dz_auc2 = _validataion_per_dz(edges[srt_idx], G_k)

for key, val in d_dz_auc2.items():
	if val > 0.7:
		print key, val
