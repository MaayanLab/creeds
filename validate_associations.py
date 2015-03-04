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
os.chdir('D:\Zichen_Projects\microtask_GEO\DisGeNET')
# with open ('curated_gene_disease_associations.txt') as f:
# with open ('all_gene_disease_associations.txt') as f:
with open ('befree_gene_disease_associations.txt') as f:
	next(f)
	for line in f:
		sl = line.strip().split('\t')
		gene, umls = sl[1], sl[3].split(':')[1]
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
		else:
			fp += 1
		tps.append(tp)
		fps.append(fp)
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
	

def plot_roc(edgelist_fn, d_gs, ax, ls='-', label=None):
	scores = np.array(map(float, file2list(edgelist_fn,2)))
	edges = np.array( zip(file2list(edgelist_fn,0), file2list(edgelist_fn, 1)) )

	fps1, tps1 = _validation2(edges, d_gs)
	srt_idx = abs(scores).argsort()[::-1]
	ranks1 = rankdata(scores, method='min')[::-1]
	FPR1, TPR1 = ranks1/max(ranks1), tps1/float(tps1[-1])
	roc_auc1 = auc(FPR1, TPR1)
	ax.plot(FPR1, TPR1, color=COLORS10[0], label='signed Jaccard' +' (AUC = %0.3f)' % roc_auc1, linewidth=2, ls=ls)

	fps2, tps2 = _validation2(edges[srt_idx], d_gs)
	abs_scores_sorted = abs(scores[srt_idx])
	ranks2 = rankdata(abs_scores_sorted, method='min')[::-1]
	FPR2, TPR2 = ranks2/max(ranks2), tps2/float(tps2[-1])
	roc_auc2 = auc(FPR2, TPR2)
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


os.chdir('D:\Zichen_Projects\microtask_GEO')
fig = plt.figure()
ax = fig.add_subplot(111)

PPIN_fns = [
'D:\Zichen_Projects\\bias_correction\PPIN\PPIN_to_months_low.p',
'D:\Zichen_Projects\\bias_correction\PPIN\PPIN_to_months.p',
]

# plot_roc('gene_umls_edgelist_signed_jaccard.txt', d_umls_gene , ax)

plot_roc('gene_gene_edgelist_signed_jaccard.txt', PPIN_fns[0] , ax, ls='-' )
plot_roc('gene_gene_edgelist_signed_jaccard.txt', PPIN_fns[1] , ax, ls='--' )

fig.tight_layout()
plt.show()
