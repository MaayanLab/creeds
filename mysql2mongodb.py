## gather cleaned data from the mysql db to mongo db
## created on 6/29/2015

import os, sys
import MySQLdb
import MySQLdb.cursors
import pymongo
from pymongo import MongoClient
sys.path.append('/Users/zichen/Documents/bitbucket/maayanlab_utils')
from fileIO import mysqlTable2dict


client = MongoClient('mongodb://127.0.0.1:27017/')
# client = MongoClient()
db = client['microtask_signatures']
coll = db['signatures']

mysql = MySQLdb.connect(host = "localhost", user = "root", passwd = "", 
	db = "maaya0_crowdsourcing")## ,cursorclass=MySQLdb.cursors.DictCursor)
cur = mysql.cursor()


d_gds_gse = mysqlTable2dict('maaya0_crowdsourcing', 'gds_gse', 0, 1)
d_gse_platform = mysqlTable2dict('maaya0_crowdsourcing', 'geo_platforms', 0, 1)
ORGANISMS = {
	'human': 'human', 'homo sapiens': 'human',
	'mouse': 'mouse', 'mus musculus': 'mouse',
	'rat': 'rat', 'rattus norvegicus': 'rat'
}

def sanitize_row(row):
	## sanitize rows in geo2enrichr* tables
	## this is generic, table specific field should be 
	## manually added to doc
	doc = {}
	geo_id, ctrls, perts, gene, pert_type, organism, cell = [item.strip() for item in row[0:7]]
	geo_id = geo_id.strip().split('[')[0]
	curator = row[-3]

	if geo_id in d_gds_gse:
		geo_id = d_gds_gse[geo_id]
	if geo_id in d_gse_platform:
		platform = d_gse_platform[geo_id]
	else:
		try:
			platform = geo_id2platform(geo_id)
		except:
			platform = 'error'
			print uid, row
			pass
	
	organism = organism.lower()
	if platform != 'error' and organism in ORGANISMS:
		organism = ORGANISMS[organism]
		ctrls = filter(None,ctrls.split(','))
		perts = filter(None,perts.split(','))
		if len(ctrls) == 1 or len(perts) == 1: 
			if ';' in ctrls[0] or ';' in perts[0]:
				delim = ';'
			else:
				delim = ' '
			ctrls = [s.strip() for s in ctrls[0].split(delim)]
			perts = [s.strip() for s in perts[0].split(delim)]

		if len(ctrls) != 1 and len(perts) != 1: # should have replicates
			if len(set(ctrls) & set(perts)) == 0: # GSMs are not allowed to be in both ctrls and perts
				doc['geo_id'] = geo_id
				doc['ctrl_ids'] = ctrls
				doc['pert_ids'] = perts
				doc['platform'] = platform
				doc['organism'] = organism
				doc['curator'] = curator
				doc['cell_type'] = cell
	return doc

def add_gene_entries():
	## move single gene pert entries
	d_uid_hssymbol = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_genes', 0, 1)
	d_uid_mmsymbol = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_genes', 0, 2)
	cur.execute("""SELECT * FROM geo2enrichr""")
	for row in cur:
		uid = row[-1]
		if uid in d_uid_hssymbol:
			hs_symbol = d_uid_hssymbol[uid]
			mm_symbol = d_uid_mmsymbol[uid]
			if hs_symbol == 'NULL': hs_symbol = None
			if mm_symbol == 'NULL': mm_symbol = None
			if hs_symbol is None and mm_symbol is None:
				pass
			else:
				doc = sanitize_row(row)
				if doc != {}:
					prefix_id = 'gene:' + str(uid)

					## table specific fields
					pert_type = row[4]
					doc['id'] = prefix_id
					doc['pert_type'] = pert_type
					doc['hs_gene_symbol'] = hs_symbol
					doc['mm_gene_symbol'] = mm_symbol
					# print coll.find_one({'id': prefix_id})				
					# print doc
					result = coll.insert_one(doc)

	return


def add_dz_entries():
	## move dz entries
	d_uid_dzname = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_dzs', 0, 1)
	d_uid_doid = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_dzs', 0, 2)
	d_uid_umls = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_dzs', 0, 3)
	cur.execute("""SELECT * FROM geo2enrichr_dz""")
	for row in cur:
		uid = row[-1]
		if uid in d_uid_dzname:
			doc = sanitize_row(row)
			if doc != {}:
				prefix_id = 'dz:' + str(uid)
				dzname = d_uid_dzname[uid].strip('"')
				doid = d_uid_doid[uid]
				umls = d_uid_umls[uid]
				doc['id'] = prefix_id
				doc['do_id'] = doid
				doc['umls_cui'] = umls
				doc['disease_name'] = dzname
				result = coll.insert_one(doc)
	return

def add_drug_entries():
	## move drug entries
	d_uid_drugname = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_drugs', 0, 1)
	d_uid_drugname_ori = mysqlTable2dict('maaya0_crowdsourcing', 'geo2enrichr_drug', -1, 3)
	d_uid_dbid = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_drugs', 0, 2)
	d_uid_pcid = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_drugs', 0, 3)
	d_uid_smiles = mysqlTable2dict('maaya0_crowdsourcing', 'drug_smiles', 1, 0)

	cur.execute("""SELECT * FROM geo2enrichr_drug""")
	for row in cur:
		uid = row[-1]
		if uid in d_uid_drugname:
			doc = sanitize_row(row)
			if doc != {}:
				prefix_id = 'drug:' + str(uid)
				drugname = d_uid_drugname[uid]
				if drugname is None:
					drugname = d_uid_drugname_ori[uid]

				dbid = d_uid_dbid[uid]
				pcid = d_uid_pcid[uid]
				smiles = d_uid_smiles[uid]

				doc['id'] = prefix_id
				doc['drugbank_id'] = dbid
				doc['pubchem_cid'] = pcid
				doc['drug_name'] = drugname
				doc['smiles'] = smiles
				result = coll.insert_one(doc)
	return


coll.drop()
print coll.count()

coll.create_index('id', unique=True)

add_gene_entries()
add_dz_entries()
add_drug_entries()


print coll.count()

mysql.close()

