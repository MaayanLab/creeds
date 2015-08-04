## ORMs and utils for computing similarity and etc.

import os, sys, json
import time
import numpy as np
from pymongo import MongoClient

# client = MongoClient('mongodb://127.0.0.1:27017/')
client = MongoClient('mongodb://146.203.54.131:27017/')
db = client['microtask_signatures']
COLL = db['signatures']
ALL_UIDS = COLL.distinct('id')


import multiprocessing
## the following two functions are used to avoid pickling error when using multiprocessing.map
## http://stackoverflow.com/questions/3288595/multiprocessing-using-pool-map-on-a-function-defined-in-a-class
def fun(f,q_in,q_out):
    while True:
        i,x = q_in.get()
        if i is None:
            break
        q_out.put((i,f(x)))

def parmap(f, X, nprocs = multiprocessing.cpu_count()):
	q_in   = multiprocessing.Queue(1)
	q_out  = multiprocessing.Queue()
	proc = [multiprocessing.Process(target=fun,args=(f,q_in,q_out)) for _ in range(nprocs)]
	for p in proc:
		p.daemon = True
		p.start()
	sent = [q_in.put((i,x)) for i,x in enumerate(X)]
	[q_in.put((None,None)) for _ in range(nprocs)]
	res = [q_out.get() for _ in range(len(sent))]
	[p.join() for p in proc]
	return [x for i,x in sorted(res)]


def jaccard(l1, l2):
	## case insensitive jaccard
	## l1, l2 are lists of tuples
	s1, s2 = set(map(lambda x: x[0].upper(), l1)), set(map(lambda x: x[0].upper(), l2))
	up = len(s1 & s2)
	dn = len(s1 | s2)
	if dn == 0: # to prevent 0 division error
		return 0
	else:
		return float(up)/dn

def signed_jaccard(s1, s2):
	## signed jaccard index for signatures
	j1 = jaccard(s1.up_genes, s2.up_genes)
	j2 = jaccard(s1.dn_genes, s2.dn_genes)
	j3 = jaccard(s1.dn_genes, s2.up_genes)
	j4 = jaccard(s1.up_genes, s2.dn_genes)
	return (j1 + j2 - j3 - j4) / 2

def signed_jaccard_against_doc(gs, uid):
	## gs is a Signature instance, uid is an id of a doc in the mongodb
	sig = DBSignature(uid)
	if sig.has_chdir():
		sig.fill_top_genes()
		return signed_jaccard(gs, sig)
	else:
		return None

def find_name(doc):
	## find the name for a doc in the mongodb based on uid
	uid = doc['id']
	prefix = uid.split(':')[0]
	if prefix == 'gene':
		return doc['hs_gene_symbol']
	elif prefix == 'dz':
		return doc['disease_name']
	else:
		return doc['drug_name']


class Signature(object):
	def __init__(self, name=None, meta=None, up_genes=None, dn_genes=None):
		## defaults:
		if name is None: name = ''
		if meta is None: meta = {}
		if up_genes is None: up_genes = []
		if dn_genes is None: dn_genes = []

		self.name = name
		self.meta = meta
		self.up_genes = up_genes
		self.dn_genes = dn_genes

	def calc_all_scores(self, nprocs=4):
		## calcuated signed jaccard score for this signatures against 
		## all signatures in the database
		uid_scores = parmap(lambda uid: (uid, signed_jaccard_against_doc(self, uid)), ALL_UIDS, nprocs=nprocs)
		uid_scores = [(uid, score) for uid, score in uid_scores if score] # filter out None
		return dict(uid_scores)

class DBSignature(Signature):
	## signature from mongodb
	def __init__(self, uid, projection={'_id':False, 'limma':False, 'fold_changes':False}):
		## the constructor also act as a way to query mongodb using
		## the id and return desirable fields by specifying projection
		doc = COLL.find_one({'id':uid}, projection)
		name = find_name(doc)
		if 'chdir' in doc:
			chdir = doc['chdir']
			del doc['chdir']
			self.chdir = chdir
		Signature.__init__(self, name, doc)

	def has_chdir(self):
		## assert if a signature has the chdir field
		if hasattr(self, 'chdir'): return True
		else: return False

	def fill_top_genes(self, cutoff=600):
		## get top up/dn genes based on a rank cutoff
		for gene, val in zip(self.chdir['genes'], self.chdir['vals'])[:cutoff]:
			if val > 0: self.up_genes.append( (gene, val) )
			else: self.dn_genes.append( (gene, val) )

	def to_json(self, meta_only=False):
		## to export the document into json
		json_data = self.meta
		if not meta_only:
			json_data['up_genes'] = self.up_genes
			json_data['down_genes'] = self.dn_genes
		return json.dumps(json_data)

	def calc_all_scores(self, nprocs=4, cutoff=600):
		## calcuated signed jaccard score for this signatures against 
		## all signatures in the database
		self.fill_top_genes(cutoff)
		uid_scores = parmap(lambda uid: (uid, signed_jaccard_against_doc(self, uid)), ALL_UIDS, nprocs=nprocs)
		uid_scores = [(uid, score) for uid, score in uid_scores if score] # filter out None
		return dict(uid_scores)

	def get_gene_vals(self, genes, na_val=0):
		## retrieve the values of a given list of genes
		if self.has_chdir():
			vals = []
			genes_upper = map(lambda x: x.upper(), self.chdir['genes'])
			for gene in genes:
				if gene in genes_upper:
					idx = genes_upper.index(gene)
					val = self.chdir['vals'][idx]
				else:
					val = na_val
				vals.append(val)
		return vals



def get_matrix(uids, genes, na_val=0):
	## retrieve a matrix based on uids of signatures and genes
	mat = np.zeros((len(genes), len(uids)))
	projection ={'id':True, '_id':False, 'chdir': True,
		'hs_gene_symbol':True, 'disease_name':True, 'drug_name':True}

	for j, uid in enumerate(uids):
		sig = DBSignature(uid, projection=projection)
		vals = sig.get_gene_vals(genes, na_val=na_val)
		mat[:, j] = vals
	return mat


## test
# '''
# gs = DBSignature('gene:24')
# gs.get_top_genes(600)
# print gs.to_json()

# t0 = time.time()

# d_uid_scores = gs.calc_all_scores(4, 600)
# print d_uid_scores.items()[:5]

# tt = time.time()
# print(len(d_uid_scores))
# print('time passed:', tt-t0)

# mat = get_matrix(['gene:24', 'dz:114'], ['TPM3', 'TNNT1', 'MYL2', 'ATP2A1'])
# print mat

# search_string = 'tp5'
# search_dict = {
# 	"$or":[
# 		{'hs_gene_symbol' : {"$regex": search_string,"$options":"i"}},
# 		{'mm_gene_symbol' : {"$regex": search_string,"$options":"i"}},
# 		# 'disease_name' : {"$regex": search_string},
# 		# 'drug_name' : {"$regex": search_string},
# 	]
# 	}
# docs = COLL.find(search_dict)
# print docs.count()
# '''