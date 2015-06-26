## to get structures for drugs using drugbank_id and pubchem_cid
## created on 6/26/2015
import os, sys
import pybel
from pubchempy import *

sys.path.append('C:\Users\Zichen\Documents\\bitbucket\maayanlab_utils')
from fileIO import file2list, mysqlTable2dict


def load_drugbank_smiles():
	mols = pybel.readfile('sdf', 'D:\Zichen_Projects\microtask_GEO\DrugBank\\all.sdf')
	d_dbid_smiles = {}
	for mol in mols:
		dbid = mol.data['DRUGBANK_ID']
		smiles = mol.write('can').split('\t')[0]
		d_dbid_smiles[dbid] = smiles
	return d_dbid_smiles

def pcid2smiles(pcid):
	pcid = int(pcid)
	try:
		c = Compound.from_cid(pcid)
		smiles = c.canonical_smiles
	except NotFoundError:
		smiles = None

	return smiles


d_uid_dbid = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_drugs', 0, 2)
d_uid_pcid = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_drugs', 0, 3)

d_dbid_smiles = load_drugbank_smiles()


## insert drug structures to the mysql database
# print pcid2smiles(5291)
import MySQLdb
conn = MySQLdb.connect(host='localhost',user='root', passwd='',db='maaya0_crowdsourcing')
cur = conn.cursor()

for uid in d_uid_dbid:
	dbid = d_uid_dbid[uid]
	pcid = d_uid_pcid[uid]

	if pcid is None and dbid is None:
		smiles = None
	elif pcid is not None:
		smiles = pcid2smiles(pcid)
	elif dbid is not None:
		if dbid in d_dbid_smiles:
			smiles = d_dbid_smiles[dbid]
		else:
			smiles = None
	cur.execute("""INSERT INTO drug_smiles VALUES(%s, %s)""", (uid, smiles))

conn.commit()
conn.close()

## calculate fingerprint similarity
## done by the rdkit Docker image
