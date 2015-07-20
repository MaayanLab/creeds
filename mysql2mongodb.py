## gather cleaned data from the mysql db to mongo db
## created on 6/29/2015

import os, sys, json
import MySQLdb
import MySQLdb.cursors
import pymongo
from pymongo import MongoClient

sys.path.append('../maayanlab_utils')
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

def clean_genes2(genes, vals): 
	## to split /// in genes and make sure genes are unique
	## replace '.' with '_'
	cleaned = []
	cleaned_vals = []
	for gene, val in zip(genes, vals):
		gene = gene.replace('.', '_')
		if '///' in gene:
			gs = gene.split('///')
			if len(set(gs) - set(cleaned)) != 0: # has non-overlaping genes
				cleaned.append(list(set(gs) - set(cleaned))[0])
				cleaned_vals.append(val)
		else:
			if gene not in cleaned:
				cleaned.append(gene)
				cleaned_vals.append(val)
	return cleaned, cleaned_vals

def add_signature_data(doc):
	## to retrieve genes from json files given doc with id filled
	d_prefix_path = {'dz':'output/microtask_dz_jsons', 
		'drug':'output/microtask_drug_jsons', 
		'gene':'output/microtask_gene_jsons'}
	prefix, uid = doc['id'].split(':')
	filepath = d_prefix_path[prefix] + '/' + uid + '.json'
	if os.path.isfile(filepath):
		data = json.load(open(filepath, 'rb'))
		up_genes = data['up_genes']
		down_genes = data['down_genes']
		genes = up_genes.keys() + down_genes.keys()
		vals = up_genes.values() + down_genes.values()
		vals = map(float, vals) # cast to float, originally unicode
		genes, vals = clean_genes2(genes, vals) ## to split /// in genes and make sure genes are unique
		doc['chdir'] = dict(zip(genes, vals))
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
					doc = add_signature_data(doc)
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
				doc = add_signature_data(doc)
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
				doc = add_signature_data(doc)
				result = coll.insert_one(doc)
	return


# coll.drop()
# print coll.count()

# coll.create_index('id', unique=True)

# add_gene_entries()
# add_dz_entries()
# add_drug_entries()

# print coll.count()

mysql.close()

print coll.count()

## update 'chdir' field
# def get_sorted_chdir_in_doc(uid):
# 	doc = coll.find_one({'id': uid})
# 	new_chdir = None
# 	if 'chdir' in doc:
# 		new_chdir = {}
# 		chdir = doc['chdir']
# 		sorted_chdir = sorted(chdir.items(), key=lambda x: abs(x[1]), reverse=True)
# 		new_chdir['genes'] = map(lambda x: x[0], sorted_chdir)
# 		new_chdir['vals'] = map(lambda x: x[1], sorted_chdir)
# 	return new_chdir


# all_uid = coll.distinct('id')
# i = 0
# for uid in all_uid:
# 	new_chdir = get_sorted_chdir_in_doc(uid)
# 	if new_chdir is not None:
# 		coll.update_one({'id': uid}, {'$set': {'chdir': new_chdir}})
# 	i += 1
# 	if i % 200 == 0:
# 		print i, len(all_uid)


## update database with limma and fold_changes

def sort_limma_json(json_fn):
	json_data = json.load(open(json_fn, 'rb'))
	limma_dict = dict(json_data['limma'])
	limma_sorted = sorted(limma_dict.items(), key=lambda x: (x[1]), reverse=False) # small to large
	genes = map(lambda x:x[0], limma_sorted)
	vals = map(lambda x:x[1], limma_sorted)
	genes, vals = clean_genes2(genes, vals) ## to split /// in genes and make sure genes are unique
	limma ={'genes': genes, 'vals': vals}

	fold_changes_dict = dict(json_data['fold_changes'])
	fold_changes_sorted = sorted(fold_changes_dict.items(), key=lambda x: (x[1]), reverse=True) #  large to small
	genes = map(lambda x:x[0], fold_changes_sorted)
	vals = map(lambda x:x[1], fold_changes_sorted)
	genes, vals = clean_genes2(genes, vals) ## to split /// in genes and make sure genes are unique
	fold_changes ={'genes': genes, 'vals': vals}	

	return limma, fold_changes

all_uid = coll.distinct('id')
# i = 0
# d_prefix_path = {'dz':'microtask_dz_jsons_limma', 'drug':'microtask_drug_jsons_limma'}
# for uid in all_uid:
# 	prefix, id = uid.split(':')
# 	json_fn = 'output/microtask_%s_jsons_limma/%s.json' %(prefix, id)
# 	if os.path.isfile(json_fn):
# 		limma, fold_changes = sort_limma_json(json_fn)
# 		coll.update_one({'id': uid}, {'$set': {'limma': limma}})
# 		coll.update_one({'id': uid}, {'$set': {'fold_changes': fold_changes}})
# 	i += 1
# 	if i % 200 == 0:
# 		print i, len(all_uid)


## correction fold change: it was calculated as A.mean(axis=1) / B.mean(axis=1), which was wrong
import numpy as np
i = 0
d_prefix_path = {'dz':'microtask_dz_jsons_limma', 'drug':'microtask_drug_jsons_limma'}
for uid in all_uid:
	doc = coll.find_one({'id': uid}, {'fold_changes':True, '_id':False})
	if 'fold_changes' in doc:
		fold_changes = doc['fold_changes']
		fold_changes['vals'] = 1./np.array(fold_changes['vals'])
		fold_changes['vals'] = fold_changes['vals'].tolist()
		coll.update_one({'id': uid}, {'$set': {'fold_changes': fold_changes}})
	i += 1
	if i % 200 == 0:
		print i, len(all_uid)

print coll.count()
