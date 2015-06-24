## to clean drug_id in the `geo2enrichr_drug` table
## created on 6/23/2015

import os, sys
import csv
import MySQLdb


cleaned_drugs = {} # {uid: [drugbank_id, pubchem_id]}
ORGANISMS = ['human', 'mouse', 'rat', 'rattus norvegicus', 'homo sapiens', 'mus musculus']

d_dbid_pcid = {}
d_pcid_dbid = {}
d_dbid_name = {}
with open ('D:\Zichen_Projects\microtask_GEO\DrugBank\drug_links_pcid.csv') as f:
	next(f)
	for sl in csv.reader(f):
		dbid, name, pcid = sl
		if pcid != '':
			pcid = int(pcid)
			d_dbid_pcid[dbid] = pcid
			d_pcid_dbid[pcid] = dbid
		d_dbid_name[dbid] = name

print len(d_dbid_pcid), len(d_dbid_name), len(d_dbid_name)

conn = MySQLdb.connect(host='localhost',user='root', passwd='',db='maaya0_crowdsourcing')
cur = conn.cursor()

query = """SELECT * FROM geo2enrichr_drug"""
cur.execute(query)

for row in cur:
	drug_id = row[4].strip()
	organism = row[5].lower()
	uid = row[-1]
	if organism in ORGANISMS:
		meta = {'dbid':None, 'pcid': None, 'name': None}
		if drug_id.startswith('DB'):
			dbid = drug_id
			meta['dbid'] = dbid
			if dbid in d_dbid_name:
				meta['name'] = d_dbid_name[dbid]
			if dbid in d_dbid_pcid:
				meta['pcid'] = d_dbid_pcid[dbid]

		elif drug_id.startswith('CID'):
			pcid = int(drug_id.strip('CID').strip())
			meta['pcid'] = pcid
			if pcid in d_pcid_dbid:
				meta['dbid'] = d_pcid_dbid[pcid]
				meta['name'] = d_dbid_name[meta['dbid']]

		else: # pcid
			try:
				pcid = int(drug_id)
				meta['pcid'] = pcid
				if pcid in d_pcid_dbid:
					meta['dbid'] = d_pcid_dbid[pcid]
					meta['name'] = d_dbid_name[meta['dbid']]
			except ValueError:
				print drug_id
				pass
		if meta['dbid'] is not None or meta['pcid'] is not None:
			cleaned_drugs[uid] = meta

print 'Number of cleaned drug entries:', len(cleaned_drugs)

for uid in cleaned_drugs:
	meta = cleaned_drugs[uid]
	cur.execute("""INSERT INTO cleaned_drugs VALUES(%s, %s, %s, %s)""", (uid, meta['name'], meta['dbid'] , meta['pcid']))

conn.commit()
conn.close()

