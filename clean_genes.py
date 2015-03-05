## to clean gene field of entries in GEO collected from microtask1
## copy valid entries to a new sqlite db, report invalid genes
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 ## Output Type 3 (Type3) or Type 42 (TrueType)
rcParams['font.sans-serif'] = 'Arial'

import urllib2
import json
import MySQLdb
import cPickle as pickle
from itertools import combinations
from collections import Counter
from pprint import pprint

from GEOentry import *
sys.path.append('C:\Users\Zichen\Documents\\bitbucket\maayanlab_utils')
from fileIO import file2list, mysqlTable2dict, sqliteTable2dict
# sys.path.append('C:\Users\Zichen\Documents\GitHub')
# from clustergram import clustergram


GENE_SYMBOL_DB = 'gene_symbols.db'
# get dicts from sqlite db
global d_hsSynos_symbol, d_mmSynos_symbol, d_hsNames_symbol, d_mmNames_symbol, hs_symbols, mm_symbols, d_mm_hs, d_hs_mm

## synonyms to symbol
d_hsSynos_symbol = sqliteTable2dict(GENE_SYMBOL_DB, """SELECT * FROM hgnc_synonyms""", 1, 0)
sql = """SELECT mgi.mgi_symbol, mgi_synonyms.synonym FROM mgi_synonyms INNER JOIN mgi ON mgi.mgi_id=mgi_synonyms.id WHERE mgi.type='Gene'"""
d_mmSynos_symbol = sqliteTable2dict(GENE_SYMBOL_DB, sql, 1, 0)
## names to symbol
d_hsNames_symbol = sqliteTable2dict(GENE_SYMBOL_DB, """SELECT * FROM hgnc""", 2, 1 )
d_mmNames_symbol = sqliteTable2dict(GENE_SYMBOL_DB, """SELECT * FROM mgi WHERE type='Gene'""", 2, 1 )
## all valid symbols
hs_symbols = d_hsNames_symbol.values()
mm_symbols = d_mmNames_symbol.values()
## homologene
d_mm_hs = sqliteTable2dict(GENE_SYMBOL_DB, """SELECT * FROM homologene""", 2, 1)
d_hs_mm = sqliteTable2dict(GENE_SYMBOL_DB, """SELECT * FROM homologene""", 1, 2)


## data of microtask1 submissions
d_uid_gene = mysqlTable2dict('maaya0_crowdsourcing', 'geo2enrichr', -1, 3)
d_uid_sp = mysqlTable2dict('maaya0_crowdsourcing', 'geo2enrichr', -1, 5)

## clean genes
global d_name_symbol, d_syno_symbol, all_symbols
d_name_symbol = dict(d_mmNames_symbol.items() + d_hsNames_symbol.items())
d_syno_symbol = dict(d_mmSynos_symbol.items() + d_hsSynos_symbol.items())
all_symbols = set(mm_symbols) | set(hs_symbols)

def clean_genes(gene):

	if gene in all_symbols:
		symbol = gene
	else:
		if gene in d_syno_symbol:
			symbol = d_syno_symbol[gene]
		elif gene in d_name_symbol:
			symbol = d_name_symbol[gene]
		else:
			symbol = "NULL"
	return symbol

def get_homologene(symbol):
	hs_symbol, mm_symbol = "NULL", "NULL"
	if symbol in hs_symbols:
		hs_symbol = symbol
		if symbol in d_hs_mm:
			mm_symbol = d_hs_mm[symbol]
	
	elif symbol in mm_symbols:
		mm_symbol = symbol
		if symbol in d_mm_hs:
			hs_symbol = d_mm_hs[symbol]
	
	return hs_symbol, mm_symbol

d_sp = {'Homo sapiens': 'human', 'Mus musculus': 'mouse', 'Rattus norvegicus':'rat'}


c = 0 # count double NULL
# with open ('cleaned_genes.txt', 'w') as out:
# 	for uid, gene_raw in d_uid_gene.items():
# 		gene = gene_raw.strip()
# 		if '(' in gene:
# 			gene = gene.split('(')[0].strip()

# 		sp = d_uid_sp[uid]
# 		symbol = clean_genes(gene)
# 		hs_symbol, mm_symbol = get_homologene(symbol)
# 		if hs_symbol == "NULL" and mm_symbol == "NULL":
# 			c += 1
# 			print uid, gene_raw
# 		out.write(str(uid) + '\t' + gene_raw + '\t' + hs_symbol + '\t' + mm_symbol + '\n')
# print 'entries with genes cannot be mapped:' , c

## add human or mouse symbols for manually cleaned file
out = open('cleaned_genes.manual.full.txt', 'w') ## then insert this file to a table `cleaned_genes` in maaya0_crowdsourcing
with open ('cleaned_genes.manual.txt') as f:
	for line in f:
		sl = line.strip().split('\t')
		uid, gene_raw, hs_symbol, mm_symbol = sl
		if hs_symbol == "NULL" and mm_symbol != "NULL":
			hs_symbol, mm_symbol = get_homologene(mm_symbol)
		elif hs_symbol != "NULL" and mm_symbol == "NULL":
			hs_symbol, mm_symbol = get_homologene(hs_symbol)
		out.write('\t'.join([uid, hs_symbol, mm_symbol]) + '\t\n')
out.close()
