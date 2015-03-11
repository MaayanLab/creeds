## to analyze entries in GEO collected with the local geo2enrichr API
## created on 3/11/2015
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 ## Output Type 3 (Type3) or Type 42 (TrueType)
rcParams['font.sans-serif'] = 'Arial'

import urllib2, json, MySQLdb
import cPickle as pickle
from itertools import combinations
from collections import Counter
from pprint import pprint

from GEOentry import *
sys.path.append('C:\Users\Zichen\Documents\\bitbucket\maayanlab_utils')
from fileIO import file2list, mysqlTable2dict, sqliteTable2dict
sys.path.append('C:\Users\Zichen\Documents\GitHub')
from clustergram import clustergram


DATADIR = 'output/annot_jsons/'
BASE_URL = 'http://127.0.0.1:8083/g2e/full?'
HOMOLOGENE_DB = 'gene_symbols.db'

## get dicts from sqlite
d_mm_hs = sqliteTable2dict(HOMOLOGENE_DB, 'SELECT * FROM homologene', 2,1) 
d_hs_ = sqliteTable2dict(HOMOLOGENE_DB, 'SELECT * FROM hgnc', 1,2) 
lm1000 = file2list('rid_lm1000.txt',1) # the L1000 genes

gene_entries = []
os.chdir(DATADIR)

## find invalid jsons and record them
# fns = os.listdir(os.getcwd())
# c = 0
# error_log = open('no_genes_errors.log', 'w')
# for fn in fns:
# 	if fn.endswith('.json'):
# 		entry = json2entry(fn, meta_only=True)
# 		if entry.chdir == 0:
# 			line = [entry.uid, entry.status, entry.message, entry.failed_to_download] 
# 			error_log.write('\t'.join(map(str, line)) + '\n')
# 			c += 1
# print 'errors:', c
# error_log.close()


## get platforms without annot files
# platforms = file2list('no_genes_errors.log', -1)
# print len(set(platforms))
# pprint(set(platforms))
# pprint(dict(Counter(platforms)))

# print len([p for p in platforms if p.startswith('GPL')])
# entry = json2entry('1014.json', meta_only=True)
# print entry.get_url(BASE_URL)
# print entry.get_url()

## get valid entries from `valid_microtask1` table and fill in data from json files
d_uid_geoid = mysqlTable2dict('maaya0_crowdsourcing', 'valid_microtask1', 0, 1)
# for fn in fns:
# 	if fn.endswith('.json'):
# 		uid = int(fn.split('.')[0])
# 		if uid in d_uid_geoid:
# 			try:
# 				entry = json2entry(fn, meta_only=True)
# 				# if entry.chdir == 0: # not supported platforms
# 				# 	print entry.uid
# 				if entry.chdir != 0:
# 					gene_entries.append(entry)
# 			except:
# 				print fn
# pickle.dump(gene_entries, open('valid_gene_entries_meta.p', 'wb'))

gene_entries = pickle.load(open('valid_gene_entries_meta.p', 'rb'))
print 'number of valid gene entries:', len(gene_entries)

## hist for number of genes measured
# plt.hist([e.chdir for e in gene_entries])
# plt.xlabel('Genes measured', fontsize=18)
# plt.show()

## check if the platform all have the same set of genes
## anwser: NO!
# for e in gene_entries:
	# if e.chdir < 5000:
	# 	print e.platform, e.chdir
	# if e.platform == 'GPL1261':
	# 	print e.chdir


## find the minimum set of genes that are shared in all entries
nums_homologenes = []
nums_lmgenes = []

print len(set(lm1000) & set(d_hs_))
print set(lm1000) - set(d_hs_)

for e in gene_entries:
	fn = str(e.uid)+'.json'
	genes = json2genes(fn)
	if e.organism not in ['human', 'Homo sapiens']:
		# genes = [d_mm_hs[gene] for gene in genes if gene in d_mm_hs]
		genes_humanize = []
		for gene in genes:
			if gene in d_mm_hs:
				genes_humanize.append(d_mm_hs[gene])
			else:
				genes_humanize.append(gene.upper())
		genes = genes_humanize
# 	homologenes = max( len(set(genes) & set(d_mm_hs.keys())), len(set(genes) & set(d_mm_hs.values())) ) #
# 	nums_homologenes.append(homologenes)
	lmgenes = len(set(genes) & set(lm1000))
	nums_lmgenes.append(lmgenes)

print len(d_mm_hs)
print len(set(d_mm_hs.values()) & set(lm1000))

print set(lm1000) - set(d_mm_hs.values())


print max(nums_lmgenes)
plt.hist(nums_lmgenes, bins=50)
plt.xlabel('Number of LM1000', fontsize=18)
plt.show()
