## script to parse file from diseaes-ontology 
## 'HumanDO.obo'
## created on 3/20/2015
## used for cleaning disease_id

import os, sys
sys.path.append('C:\\Users\\Zichen\\Documents\\bitbucket\\maayanlab_utils')
import obo_parser


DO_FILE = 'D:\Zichen_Projects\microtask_GEO\DO\HumanDO.obo'
ontology_id_terms, name2id, G = obo_parser.obo2network(DO_FILE)
print len(ontology_id_terms), len(name2id), G.number_of_edges(), G.number_of_nodes()
# print ontology_id_terms['DOID:0001816'].xrefs
# print ontology_id_terms['DOID:0001816'].synonym
# print ontology_id_terms['DOID:0001816'].name
# print ontology_id_terms['DOID:0001816'].parents
# print ontology_id_terms['DOID:0001816'].children

d_doid_umls_ids = {}
d_umls_doid = {}

for doid, term in ontology_id_terms.items():
	umls_ids = [xref for xref in term.xrefs if xref.startswith('UMLS_CUI:')]
	d_doid_umls_ids[doid] = umls_ids
	for umls_id in umls_ids:
		d_umls_doid[umls_id] = doid

# output the dictionary
# with open ('umls_id2_DOID.info', 'w') as out:
# 	for umls_id, doid in d_umls_doid.items():
# 		out.write(umls_id+'\t'+doid+'\n')

## map between DOID and UMLS_CUI
from fileIO import mysqlTable2dict
d_uid_dzid = mysqlTable2dict('maaya0_crowdsourcing', 'geo2enrichr_dz', -1, 4)
cleaned_dzids = {}
for uid, dzid in d_uid_dzid.items():
	if ':' in dzid:
		dzid = ':'.join([s.strip() for s in dzid.split(':') if s])
	if dzid.startswith('DOID:'):
		doid = dzid
		umls_ids = d_doid_umls_ids[doid]
	elif dzid.startswith('C'):
		umls_ids = ['UMLS_CUI:' + dzid]
		if 'UMLS_CUI:'+dzid in d_umls_doid: doid = d_umls_doid['UMLS_CUI:'+dzid]
		else: doid = None
	elif dzid.startswith('UMLS_CUI'):
		umls_ids = [dzid]
		if dzid in d_umls_doid: doid = d_umls_doid[dzid]
		else: doid = None
	else:
		doid = 'DOID:' + dzid
		umls_ids = d_doid_umls_ids[doid]
	cleaned_dzids[uid] = [doid, umls_ids]	

print cleaned_dzids[650]

## insert disease id to the mysql database
import MySQLdb
conn = MySQLdb.connect(host='localhost',user='root', passwd='',db='maaya0_crowdsourcing')
cur = conn.cursor()

for uid in cleaned_dzids:
	doid, umls_ids = cleaned_dzids[uid]
	if len(umls_ids) == 0:
		cur.execute("""INSERT INTO cleaned_dzs VALUES(%s, %s, %s)""", (uid, doid, None))
	else:	
		for umls_id in umls_ids:
			umls_id = umls_id.split(':')[1]
			cur.execute("""INSERT INTO cleaned_dzs VALUES(%s, %s, %s)""", (uid, doid, umls_id))

conn.commit()
conn.close()
