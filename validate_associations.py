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

import json
from GEOentry import *
sys.path.append('C:\Users\Zichen\Documents\\bitbucket\maayanlab_utils')
from fileIO import file2list
from plots import enlarge_tick_fontsize, COLORS10
from homologene import *

## parse do_rif

d_geneID_geneSymbol = id2symbol('9606')

d_umls_gene = {} # to store dz-gene associations from do_rif.human

# os.chdir('D:\Zichen_Projects\microtask_GEO\do_rif')
# with open ('do_rif.human.txt') as f:
# 	for line in f:
# 		sl = line.strip().split('\t')
# 		geneID, umls = sl[0], sl[3]
# 		if geneID in d_geneID_geneSymbol:
# 			gene = d_geneID_geneSymbol[geneID]
# 			if umls not in d_umls_gene:
# 				d_umls_gene[umls] = [gene]
# 			else:
# 				if gene not in d_umls_gene[umls]:
# 					d_umls_gene[umls].append(gene)
# os.chdir('D:\Zichen_Projects\microtask_GEO\DisGeNET')
# with open ('curated_gene_disease_associations.txt') as f:
# # with open ('all_gene_disease_associations.txt') as f:
# # with open ('befree_gene_disease_associations.txt') as f:
# 	next(f)
# 	for line in f:
# 		sl = line.strip().split('\t')
# 		gene, umls = sl[1], sl[3].split(':')[1]
# 		if umls not in d_umls_gene:
# 			d_umls_gene[umls] = [gene]
# 		else:
# 			if gene not in d_umls_gene[umls]:
# 				d_umls_gene[umls].append(gene)

os.chdir('D:\Zichen_Projects\microtask_GEO\DISEASES')
with open ('human_disease_knowledge_filtered.tsv') as f:

	for line in f:
		sl = line.strip().split('\t')
		gene, umls = sl[1], sl[2]
		if umls not in d_umls_gene:
			d_umls_gene[umls] = [gene]
		else:
			if gene not in d_umls_gene[umls]:
				d_umls_gene[umls].append(gene)


print len(d_umls_gene)


def _validation(edges, d_gs):
	fps = []
	tps = []
	fp = 0
	tp = 0
	for gene, umls in edges:
		if umls in d_umls_gene:
			if gene in d_umls_gene[umls]:
				tp += 1
			else:
				fp += 1
			tps.append(tp)
			fps.append(fp)
		else:
			# fp += 1
			continue
	max_tp = float(tps[-1])
	max_fp = float(fps[-1])
	print 'max_tp=', max_tp
	return np.array(fps), np.array(tps)

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
	

def plot_roc(edgelist_fn, d1, d2, d_gs, ax, ls='-', label=None):
	## d1: d_uid_name for col1
	## d2: d_uid_name for col2

	scores_ = file2list(edgelist_fn,2)
	nodes1 = [d1[int(s)] for s in file2list(edgelist_fn,0)]
	nodes2 = [d2[int(s)] for s in file2list(edgelist_fn,1)]
	
	edges = []
	scores = []
	for n1, n2, score in zip(nodes1, nodes2, scores_):
		if n1 != None and n2 != None:
			edges.append( (n1, n2) )
			scores.append(float(score))
	edges = np.array(edges)
	scores = np.array(scores)

	fps1, tps1 = _validation(edges, d_gs)
	srt_idx = abs(scores).argsort()[::-1]
	# ranks1 = rankdata(scores, method='min')[::-1]
	# FPR1, TPR1 = ranks1/max(ranks1), tps1/float(tps1[-1])
	FPR1, TPR1 = fps1/float(fps1[-1]), tps1/float(tps1[-1])
	roc_auc1 = auc(FPR1, TPR1)
	# roc_auc1 = 0
	ax.plot(FPR1, TPR1, color=COLORS10[0], label='signed Jaccard' +' (AUC = %0.3f)' % roc_auc1, linewidth=2, ls=ls)

	fps2, tps2 = _validation(edges[srt_idx], d_gs)
	abs_scores_sorted = abs(scores[srt_idx])
	ranks2 = rankdata(abs_scores_sorted, method='min')[::-1]
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
	ax.legend(loc='lower right',prop={'size':18})



## load dz id
d_uid_doid = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_dzs', 0, 1)
d_uid_umls = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_dzs', 0, 2)

d_uid_gene = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_genes', 0, 1)


os.chdir('C:\\Users\Zichen\Documents\\bitbucket\\microtask_GEO\\output\\analyses_output')

PPIN_fns = [
'D:\Zichen_Projects\\bias_correction\PPIN\PPIN_to_months_low.p',
'D:\Zichen_Projects\\bias_correction\PPIN\PPIN_to_months.p',
]

fig = plt.figure()
ax = fig.add_subplot(111)

# plot_roc('gene_umls_edgelist_signed_jaccard.txt', d_umls_gene , ax)

# plot_roc('gene_gene_edgelist_signed_jaccard.txt', PPIN_fns[0] , ax, ls='-' )
# plot_roc('gene_gene_edgelist_signed_jaccard.txt', PPIN_fns[1] , ax, ls='--' )

plot_roc('gene_dz_edgelist_signed_jaccard.txt', d_uid_gene, d_uid_doid, d_umls_gene, ax)
fig.tight_layout()
plt.show()
