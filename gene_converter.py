## to convert gene symbols to Entrez gene ID
import sqlite3

## this database is created based homologene, and data from HGNC and MGI
HOMOLOGENE_DB = 'static/data/gene_symbols.db'

def sqliteTable2dict(conn, query, key_idx, val_idx):
	d = {}
	cur = conn.cursor()
	cur.execute(query)
	for row in cur:
		key, val = row[key_idx], row[val_idx]
		d[key] = val
	return d

def load_gene_symbol_dict():
	conn = sqlite3.connect(HOMOLOGENE_DB)
	d_hs = sqliteTable2dict(conn, 'SELECT * FROM hgnc', 1,3) # all human symbols 
	d_mm = sqliteTable2dict(conn, """SELECT * FROM mgi WHERE type='Gene'""", 1,4) # all mouse symbols

	GENE_SYMBOLS = dict(d_hs.items() + d_mm.items())
	conn.close()
	return GENE_SYMBOLS
