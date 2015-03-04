## make MDS plot
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn import manifold

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 ## Output Type 3 (Type3) or Type 42 (TrueType)
rcParams['font.sans-serif'] = 'Arial'

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

d_uid_gene = mysqlTable2dict('maaya0_crowdsourcing', 'geo2enrichr', -1, 3)
gene_entries = []
os.chdir('D:\Zichen_Projects\microtask_GEO\jsons')
fns = os.listdir(os.getcwd())
for fn in fns:
	entry = json2entry(fn)
	if len(entry.chdir) >5000:
		entry.get_lists_cutoff(500)
		gene_entries.append(entry)

os.chdir('D:\Zichen_Projects\microtask_GEO')
mds = manifold.MDS(n_components=2, dissimilarity="precomputed")
nmds = manifold.MDS(metric=False, n_components=2, dissimilarity="precomputed", n_init=1)
dissimilarity_matrix = 1 - abs(np.loadtxt('genes_signed_Jaccard_matrix_n1663x1663.txt'))
# dissimilarity_matrix = 1- abs(np.loadtxt('dzs_signed_Jaccard_matrix_n171x171.txt'))
mds.fit(dissimilarity_matrix)
pos = mds.embedding_

npos = nmds.fit_transform(dissimilarity_matrix, init=pos)


print pos.shape, npos.shape

with open ('genes_1663_MDS_coords.tsv', 'w') as out:
	for (x, y), e in zip(npos, gene_entries):
		out.write( '\t'.join(map(str,  [x,y, e.uid, e.geo_id, d_uid_gene[e.uid].strip(), e.pert_type, e.cell, e.curator, e.platform] ) ) +'\n')

# fig = plt.figure(figsize=(10,6))
# ax1 = fig.add_subplot(121)
# ax2 = fig.add_subplot(122)

# ax1.scatter(pos[:,0], pos[:,1])
# ax2.scatter(npos[:,0], npos[:,1])
# plt.show()
