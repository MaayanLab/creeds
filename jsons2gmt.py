## to read gene signatures from json files and write to gmt
## created on 6/24/2015

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
from gene_convertion import *

sys.path.append('C:\Users\Zichen\Documents\\bitbucket\maayanlab_utils')
from fileIO import file2list, mysqlTable2dict, sqliteTable2dict, read_gmt
sys.path.append('C:\Users\Zichen\Documents\GitHub')
from clustergram import clustergram



def get_valid_entries(datadir, outfn):
	"""
	to read meta data from all the json file in a datadir
	and return the meta data dict key by uid
	"""
	# cwd = os.getcwd() # remember initial dir
	os.chdir(datadir)
	if outfn in os.listdir(os.getcwd()):
		valid_entries = pickle.load(open(outfn, 'rb'))
	else:
		valid_entries = {}
		fns = os.listdir(os.getcwd())
		for fn in fns:
			if fn.endswith('.json'):
				uid = int(fn.split('.')[0])
				try:
					entry = json2entry(fn, meta_only=True)
					if entry.chdir > 5000: # at least 5000 genes measured
						valid_entries[uid] = entry
				except:
					print fn
					pass
	pickle.dump(valid_entries, open(outfn, 'wb'))
	print 'Number of valid entries', len(valid_entries)
	return valid_entries

def make_gmt(datadir, meta_fn, gmt_fn):

	valid_entries = get_valid_entries(datadir, meta_fn)

	d_uid_umls_id = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_dzs', 0, 2)
	CUTOFF = 600
	i = 0
	with open(gmt_fn, 'w') as out:
		for uid, e in valid_entries.items():
			i += 1
			fn = str(e.uid)+'.json'
			entry = json2entry(fn, meta_only=False) # full entry
			geneset = entry.to_json_geneset2(CUTOFF)
			# geneset['desc'] = d_uid_umls_id[e.uid]
			geneset['desc'] = str(uid)

			try:
				up_genes = []
				dn_genes = []
				for gene, val in zip(geneset['genes'], geneset['vals']):
					if val > 0:
						up_genes.append(gene)
					else:
						dn_genes.append(gene)
				out.write(geneset['term'] +'|up' + '\t' + geneset['desc'] + '\t' + '\t'.join(up_genes) + '\n')
				out.write(geneset['term'] +'|dn' + '\t' + geneset['desc'] + '\t' + '\t'.join(dn_genes) + '\n')

			except:
				pass
			if i % 25 == 0:
				print i, e.gene, geneset['term']
	return


# make_gmt('output/microtask_gene_jsons', 'valid_gene_entries.pkl', 'crowdsourced_single_gene_pert_top600.gmt')
# make_gmt('output/microtask_dz_jsons', 'valid_dz_entries.pkl', 'crowdsourced_diseases_top600.gmt')
# make_gmt('output/microtask_drug_jsons', 'valid_drug_entries.pkl', 'crowdsourced_drugs_top600.gmt')

