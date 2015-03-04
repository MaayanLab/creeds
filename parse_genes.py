## parse data downloaded from HGNC and MGI
## to create a database of gene symbols and their synonyms
# HGNC txt file downloaded from http://tinyurl.com/o8mmfrx

import os, sys
import sqlite3

WORKDIR = "D:\Zichen_Projects\microtask_GEO\HGNC_and_MGI"

conn = sqlite3.connect('gene_symbols.db')
cursor = conn.cursor()
# Create tables with schema defined in the sql script
query = open('gene_symbols.sql', 'r').read()
cursor.executescript(query)
conn.close()

conn = sqlite3.connect('gene_symbols.schema.db')
cursor = conn.cursor()
## list all tables
cursor.execute('''SELECT name FROM sqlite_master WHERE type='table';''')
print cursor.fetchall()

## parse data
os.chdir(WORKDIR)
hgnc_entries = []
hgnc_synonyms = []
with open ('HGNC.txt') as f:
	next(f)
	for line in f:
		sl = line.split('\t')
		status = sl[3]
		if status == 'Approved':
			symbol = sl[1]
			if sl[8] != '\n':	entrez_id = int(sl[8])
			else: entrez_id = "NULL"

			hgnc_entry = ( sl[0], symbol, sl[2], entrez_id )
			hgnc_entries.append( hgnc_entry )

			previous_symbols = [s.strip() for s in sl[4].split(',')]
			previous_names = [s for s in filter(None, sl[5].split('"')) if s != ', ']
			synonyms = [ s.strip() for s in sl[6].split(',') ]
			name_synonyms = [s for s in filter(None, sl[7].split('"')) if s != ', ']
			all_synonyms = synonyms + name_synonyms + previous_symbols + previous_names
			all_synonyms = filter(None, all_synonyms)
			for synonym in all_synonyms:
				hgnc_synonyms.append( (symbol, synonym) )

mgi_entries = []
mgi_synonyms = []
with open ('MGI_EntrezGene.rpt') as f:
	for line in f:
		sl = line.split('\t')
		status = sl[2]
		if status != 'W':
			symbol = sl[1]
			if sl[8] != '': entrez_id = int(sl[8])
			else: entrez_id = "NULL"
			mgi_entry = (sl[0], symbol, sl[3], sl[6], entrez_id )
			mgi_entries.append( mgi_entry )
			
			synonyms = filter(None, sl[9].split('|'))
			for synonym in synonyms:
				mgi_synonyms.append( (symbol, synonym) )

d_id_symbols = {}
with open ('HOM_MouseHumanSequence.rpt') as f:
	next(f)
	for line in f:
		sl = line.split('\t')
		id = sl[0]
		symbol = sl[3]
		taxid = sl[2]
		if id not in d_id_symbols:
			d_id_symbols[id] = {'9606':None, '10090': None}
		d_id_symbols[id][taxid] = symbol

homologene_entries = []
for id in d_id_symbols:
	if d_id_symbols[id]['9606'] != None and d_id_symbols[id]['10090'] != None:
		homologene_entries.append( (id, d_id_symbols[id]['9606'], d_id_symbols[id]['10090']) )


## insert data into tables
print map(len, [hgnc_entries, hgnc_synonyms, mgi_entries, mgi_synonyms, homologene_entries])
cursor.executemany('INSERT INTO hgnc VALUES (?,?,?,?)', hgnc_entries)
cursor.executemany('INSERT INTO hgnc_synonyms VALUES (?,?)', hgnc_synonyms)
cursor.executemany('INSERT INTO mgi VALUES (?,?,?,?,?)', mgi_entries)
cursor.executemany('INSERT INTO mgi_synonyms VALUES (?,?)', mgi_synonyms)
cursor.executemany('INSERT INTO homologene VALUES (?,?,?)', homologene_entries)
conn.commit()
conn.close()
