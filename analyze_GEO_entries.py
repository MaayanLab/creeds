## to analyze entries in GEO collected
## created on 2/13/2015
import os, sys
import numpy as np
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
from collections import Counter
from pprint import pprint

from GEOentry import *
sys.path.append('C:\Users\Zichen\Documents\\bitbucket\maayanlab_utils')
from fileIO import file2list, mysqlTable2dict
sys.path.append('C:\Users\Zichen\Documents\GitHub')
from clustergram import clustergram

def list2file(l, fn):
	with open (fn, 'w') as out:
		out.write('\n'.join(l) + '\n')
	return

def index_matrix(l, func):
	## make index matrix
	mat = np.ones((len(l), len(l)))
	for i, j in combinations(range(len(l)), 2):
		mat[i, j] = func(l[i], l[j])
		mat[j, i] = func(l[j], l[i])
		if i % 50 == 0 and i != 0:
			print i, len(l)
	return mat	

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


d_uid_gene = mysqlTable2dict('maaya0_crowdsourcing', 'geo2enrichr', -1, 3)

os.chdir('D:\Zichen_Projects\microtask_GEO')
d_gse_platform = dict(zip(file2list('GSE_platform.txt',0), file2list('GSE_platform.txt', 1) ))

dz_entries = []
os.chdir('D:\Zichen_Projects\microtask_GEO\jsons_dz')
fns = os.listdir(os.getcwd())
for fn in fns:
	entry = json2entry(fn)
	if len(entry.chdir) >5000:
		entry.get_lists_cutoff(500)
		dz_entries.append(entry)

gene_entries = []
os.chdir('D:\Zichen_Projects\microtask_GEO\jsons')
fns = os.listdir(os.getcwd())
for fn in fns:
	entry = json2entry(fn)
	if len(entry.chdir) >5000:
		entry.get_lists_cutoff(500)
		gene_entries.append(entry)


os.chdir('D:\Zichen_Projects\microtask_GEO')
## output gene_entries and dz_entries
pickle.dump(gene_entries, open('gene_entries_500up500dn.p','wb'))
pickle.dump(dz_entries, open('dz_entries_500up500dn.p','wb'))

# compute matrix of Signed Jaccard
## self similarity matrix
# mat_j_gene = index_matrix(gene_entries, signed_jaccard)
# np.savetxt('genes_signed_Jaccard_matrix_n%sx%s.txt'%(len(gene_entries), len(gene_entries)), mat_j_gene, delimiter='\t')		
# del mat_j_gene

# mat_j_dz = index_matrix(dz_entries, signed_jaccard)
# np.savetxt('dzs_signed_Jaccard_matrix_n%sx%s.txt'%(len(dz_entries), len(dz_entries)), mat_j_dz, delimiter='\t')		
# del mat_j_dz

## cross gene pert with diseases
# signed_jaccard_matrix = np.ones((len(gene_entries), len(dz_entries)))
# print signed_jaccard_matrix.shape

# for i, gene_entry in enumerate(gene_entries):
# 	for j, dz_entry in enumerate(dz_entries):
# 		score = signed_jaccard(gene_entry, dz_entry)
# 		signed_jaccard_matrix[i,j] = score
# np.savetxt('signed_Jaccard_matrix_n%sx%s.txt'%(len(gene_entries), len(dz_entries)), signed_jaccard_matrix, delimiter='\t')		

## output labels
# gene_names = [d_uid_gene[e.uid].strip() for e in gene_entries]
# list2file([d_uid_gene[e.uid].strip() for e in gene_entries], 'gene_names.txt')
# list2file([e.cell for e in gene_entries], 'gene_cells.txt')
# list2file([e.curator for e in gene_entries], 'gene_curators.txt')
# list2file([e.platform for e in gene_entries], 'gene_platforms.txt')
# list2file(['|'.join([e.gene, e.pert_type, e.cell, e.organism, e.curator, e.uid ]) for e in gene_entries], 'gene_metas.txt')
# list2file(['|'.join([d_uid_gene[e.uid].strip(), e.pert_type, str(e.uid) ]) for e in gene_entries], 'gene_metas_short.txt')

# list2file([e.gene for e in dz_entries], 'dz_names.txt')
# list2file([e.cell for e in dz_entries], 'dz_cells.txt')
# list2file([e.curator for e in dz_entries], 'dz_curators.txt')
# list2file([e.platform for e in dz_entries], 'dz_platforms.txt')
# list2file(['|'.join([e.gene, e.pert_type, e.cell, e.organism, e.curator, e.uid ]) for e in dz_entries], 'dz_metas.txt')

# d_curator = {
# 	'Andrew': 'Andrew',
# 	'Kevin': 'Kevin',
# 	'Cadimo': 'Cadimo',
# 	'gszeto':'gszeto',
# 	'FabioAmaral':'FabioAmaral',
# 	'HollyWoodland':'HollyWoodland',
# 	'OliFucMuc':'OliFucMuc',
# 	'Bhanuprakash':'Bhanuprakash',
# 	'cms':'cms',
# 	'shvetank':'shvetank'
# }
# gene_curators = []
# for e in gene_entries:
# 	if e.curator in d_curator:
# 		gene_curators.append(e.curator)
# 	else:
# 		gene_curators.append('others')

# signed_jaccard_matrix = np.loadtxt('signed_Jaccard_matrix_n1663x171.txt')
# clustergram(data=signed_jaccard_matrix, row_labels=[e.platform for e in gene_entries], col_labels=[e.platform for e in dz_entries],
# 			row_groups=[e.platform for e in gene_entries], col_groups=[e.platform for e in dz_entries], cluster=True,
# 			row_linkage='average', col_linkage='average', 
# 			row_pdist='euclidean', col_pdist='euclidean',
# 			standardize=3, log=False, colormap='redbluecmap',
# 			display_range=[-.05,.05], figsize=8, figname=None, colorkey='Signed Jaccard Index')

## sort matrix and output edgelist format txt files

# edges_vals = adjMatrx2edgeList(signed_jaccard_matrix, gene_names, [e.pert_type for e in dz_entries])
# with open ('gene_umls_edgelist_signed_jaccard.txt','w') as out:
# 	for edge, val in edges_vals:
# 		out.write(edge + '\t' + str(val) + '\n')

# mat_j_gene = np.loadtxt('genes_signed_Jaccard_matrix_n1663x1663.txt')
# edges_vals = adjMatrx2edgeList(mat_j_gene, gene_names, gene_names)
# with open ('gene_gene_edgelist_signed_jaccard.txt','w') as out:
# 	for edge, val in edges_vals:
# 		out.write(edge + '\t' + str(val) + '\n')

# mat_j_dz = np.loadtxt('dzs_signed_Jaccard_matrix_n171x171.txt')
# edges_vals = adjMatrx2edgeList(mat_j_dz, [e.gene for e in dz_entries], [e.gene for e in dz_entries])
# with open ('dz_dz_edgelist_signed_jaccard.txt','w') as out:
# 	for edge, val in edges_vals:
# 		out.write(edge + '\t' + str(val) + '\n')




'''
fig = plt.figure(figsize=(16,8))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.set_title('gene perturbations', fontsize=20)
# ax1.hist([len(e.chdir) for e in gene_entries], bins=20)
ax1.scatter([len(e.chdir) for e in gene_entries], [e.conversion_pct for e in gene_entries], alpha=0.5)
ax1.set_xlabel('Number of genes measured', fontsize=16)
ax1.set_ylabel('conversion_pct', fontsize=16)

ax2.set_title('Disease signatures', fontsize=20)
# ax2.hist([len(e.chdir) for e in dz_entries], bins=20)
ax2.scatter([len(e.chdir) for e in dz_entries], [e.conversion_pct for e in dz_entries], alpha=0.5)
ax2.set_xlabel('Number of genes measured', fontsize=16)
ax2.set_ylabel('conversion_pct', fontsize=16)

plt.show()
'''

# conn = MySQLdb.connect(host='localhost',user='root', passwd='',db='maaya0_crowdsourcing')
# cur = conn.cursor()

# ## get cleaned GEO entries
# all_entries = []
# query = """SELECT * FROM geo2enrichr_cleaned"""
# cur.execute(query)
# for row in cur:
# 	geo_id, ctrls, perts, gene, pert_type, organism, cell = row[0:7]
# 	curator = row[-3]
# 	up_genes, dn_genes = row[7:9]
# 	platform = d_gse_platform[geo_id]
# 	entry = GEOentry(geo_id, ctrls.split(','), perts.split(','), gene, pert_type, platform, organism, cell, curator)
# 	entry.up_genes = np.array(up_genes.split(','))
# 	entry.dn_genes = np.array(dn_genes.split(','))
# 	all_entries.append(entry)

# print len(all_entries)
# genes = [e.gene for e in all_entries]
# perttypes = [e.pert_type for e in all_entries]
# geoids = [e.geo_id for e in all_entries]
# curators = [e.curator for e in all_entries]
# print len(set(genes)), len(set(perttypes)), len(set(geoids))

## diag matrix
# mat_j = index_matrix(all_entries) ## signed jaccard index matrix
# np.savetxt('%s_entries_signed_Jaccard_matrix.txt'%len(all_entries), mat_j, delimiter='\t')
# mat_j = np.loadtxt('%s_entries_signed_Jaccard_matrix.txt'%len(all_entries))
# list2file(genes ,'%s_entries_genes.txt'%len(all_entries))
# list2file(perttypes ,'%s_entries_perttypes.txt'%len(all_entries))
# print 'finished compute the matrix'

# clustergram(data=mat_j, row_labels=perttypes, col_labels=perttypes,
# 			row_groups=curators, col_groups=curators, cluster=True,
# 			row_linkage='average', col_linkage='average', 
# 			row_pdist='euclidean', col_pdist='euclidean',
# 			standardize=3, log=False, colormap='redbluecmap',
# 			display_range=[-.1,.1], figsize=8, figname=None, colorkey='Signed Jaccard Index')

# gene_lists = {}
# DEGs = []
# for i, e in enumerate(all_entries):
# 	gene_lists['%s_up'%i] = e.up_genes.tolist()
# 	gene_lists['%s_dn'%i] = e.dn_genes.tolist()
# 	DEGs.extend(e.up_genes.tolist())
# 	DEGs.extend(e.dn_genes.tolist())


## DEG occurrence analysis
# import matplotlib.pyplot as plt
# fig = plt.figure()
# ax = fig.add_subplot(111)
# d_gene_occ = dict(Counter(DEGs))

# with open ('DEG_occurrence.csv', 'w') as out:
# 	for gene in d_gene_occ:
# 		out.write( gene + ',' + str(d_gene_occ[gene]) + '\n')
# ax.hist(d_gene_occ.values(), log=True, bins=20)
# ax.set_ylim([0.1,1e4])
# ax.set_ylabel('# of genes')
# ax.set_xlabel('Occurrences')
# plt.show()

## gene-gene similarity matrix
# d_gene_lists = {} # gene as key and list name as values
# for list_name, genes in gene_lists.items():
# 	for gene in genes:
# 		if gene not inf d_gene_lists:
# 			d_gene_lists[gene] = [list_name]
# 		else:
# 			d_gene_lists[gene].append(list_name)

