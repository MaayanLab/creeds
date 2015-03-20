## utils used for converting gene symbols
## created on 3/20/2015

import sys
sys.path.append('C:\Users\Zichen\Documents\\bitbucket\maayanlab_utils')
from fileIO import mysqlTable2dict, sqliteTable2dict, file2list

HOMOLOGENE_DB = 'gene_symbols.db'

## get dicts from sqlite
global d_mm_hs, d_hs_
d_mm_hs = sqliteTable2dict(HOMOLOGENE_DB, 'SELECT * FROM homologene', 2,1) 
d_hs_ = sqliteTable2dict(HOMOLOGENE_DB, 'SELECT * FROM hgnc', 1,2) # all human symbols 
d_syno_hs = sqliteTable2dict(HOMOLOGENE_DB, 'SELECT * FROM hgnc_synonyms', 1,0)
lm1000 = file2list('rid_lm1000.txt',1) # the L1000 genes
global GENE_SYMBOLS
GENE_SYMBOLS = sqliteTable2dict(HOMOLOGENE_DB, """SELECT * FROM mgi WHERE type='Gene'""", 1,3).keys()
GENE_SYMBOLS += d_hs_.keys()
GENE_SYMBOLS = set(GENE_SYMBOLS)

def clean_genes(genes): ## to split /// in genes and keep only valid symbols in HGNC and MGI
	cleaned = []
	for gene in genes:
		if '///' in gene:
			for g in gene.split('///'):
				if g in GENE_SYMBOLS:
					cleaned.append(g)
		else:
			if gene in GENE_SYMBOLS:
				cleaned.append(gene)
	return cleaned

def humanize(genes): ## to convert mouse gene symbols to human's using homologene
	cleaned = clean_genes(genes)
	humanized = []
	for gene in cleaned:
		if gene in d_hs_: # human symbols already
			humanized.append(gene)
		else:
			if gene in d_mm_hs:
				humanized.append(d_mm_hs[gene])
	return humanized

