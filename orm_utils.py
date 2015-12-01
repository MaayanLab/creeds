## ORMs and utils for computing similarity and etc.

import os, sys, json
import time
import numpy as np
from pymongo import MongoClient
import requests

from gene_converter import *

## connect to mongodb
# client = MongoClient('mongodb://127.0.0.1:27017/')
client = MongoClient('mongodb://146.203.54.131:27017/')
db = client['microtask_signatures']
COLL = db['signatures']
ALL_UIDS = COLL.distinct('id')

## load gene symbol to gene ID conversion dict
GENE_SYMBOLS = load_gene_symbol_dict()

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

def load_and_fill_sig(uid):
	sig = DBSignature(uid)
	if sig.has_chdir():
		sig.fill_top_genes()
		return sig
	else:
		return None

def load_all_db_sigs(nprocs=4):
	## load all signatures with chdir from the mongodb
	d_uid_sigs = {}
	if nprocs > 1:
		sigs = parmap(load_and_fill_sig, ALL_UIDS, nprocs=nprocs)
	else:
		sigs = map(load_and_fill_sig, ALL_UIDS)
	for uid, sig in zip(ALL_UIDS, sigs):
		if sig is not None:
			d_uid_sigs[uid] = sig
	return d_uid_sigs

def find_name(doc):
	## find the name for a doc in the mongodb based on uid
	uid = doc['id']
	prefix = uid.split(':')[0]
	if prefix == 'gene':
		if doc['organism'] == 'human':
			name = doc['hs_gene_symbol']
			if name is None:
				name = doc['mm_gene_symbol']
		else:
			name = doc['mm_gene_symbol']
			if name is None:
				name = doc['hs_gene_symbol']
	elif prefix == 'dz':
		name = doc['disease_name']
	else:
		name = doc['drug_name']
	return name


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

	def calc_all_scores(self, d_uid_sigs):
		## calcuated signed jaccard score for this signatures against 
		## all signatures in d_uid_sigs
		uid_scores = []
		for uid, sig in d_uid_sigs.items():
			score = signed_jaccard(self, sig)
			uid_scores.append((uid, score))
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

	def to_dict(self, format='gmt'):
		## method to generate files for downloading
		if format == 'gmt':
			dict_data = {'name': self.name, 'id': self.meta['id']}			
		else:
			dict_data = self.meta
		dict_data['up_genes'] = self.up_genes
		dict_data['down_genes'] = self.dn_genes
		return dict_data


	def calc_all_scores(self, d_uid_sigs, cutoff=600):
		## calcuated signed jaccard score for this signatures against 
		## all signatures in the database
		self.fill_top_genes(cutoff)
		uid_scores = []
		for uid, sig in d_uid_sigs.items():
			score = signed_jaccard(self, sig)
			uid_scores.append((uid, score))
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

	def post_to_paea(self, cutoff=2000):
		## post top n genes to PAEA and return a PAEA url
		## return None if instance has no chdir
		post_url = 'http://amp.pharm.mssm.edu/Enrichr/addList'
		base_url = 'http://amp.pharm.mssm.edu/PAEA?id='
		paea_url = None
		if self.has_chdir():
			gene_list = ''
			for gene, coef in zip(self.chdir['genes'], self.chdir['vals'])[:cutoff]:
				gene_list += '%s,%s\n'% (gene, coef)
			data = {'list': gene_list, 'inputMethod': "PAEA", 'description': self.name}
			r = requests.post(post_url, files=data)
			paea_url = base_url + str(json.loads(r.text)['userListId'])
		return paea_url

	def post_to_cds2(self, cutoff=2000):
		## post top n genes to L1000CDS2 API and return a CDS2 url
		url = 'http://amp.pharm.mssm.edu/L1000CDS2/query'
		cds2_url = None
		if self.has_chdir():
			data = {
				"genes": map(lambda x: x.upper(), self.chdir['genes'][:cutoff]), 
				"vals":  self.chdir['vals'][:cutoff]
				}
			config = {"aggravate":False,"searchMethod":"CD","share":True,"combination":True,"db-version":"latest"}
			metadata = [{"key":"name","value": self.name}]
			for key, val in self.meta.items():
				if key not in ['pert_ids', 'ctrl_ids', 'curator']:
					metadata.append({"key":key, "value":val})
			payload = {"data":data,"config":config,"meta":metadata}
			headers = {'content-type':'application/json'}
			r = requests.post(url,data=json.dumps(payload),headers=headers)
			resCD = r.json()
			shareId = resCD['shareId']
			cds2_url = 'http://amp.pharm.mssm.edu/L1000CDS2/#/result/' + shareId
		return cds2_url

	def get_url(self):
		## get the url of the signature's gene, disease or drug
		url = ''
		meta = self.meta
		if meta['id'].startswith('gene:'):
			organism = meta['organism']
			if organism == 'human':
				gene_symbol = meta['hs_gene_symbol']
				if gene_symbol is None:
					gene_symbol = meta['mm_gene_symbol']
			else:
				gene_symbol = meta['mm_gene_symbol']
				if gene_symbol is None:
					gene_symbol = meta['hs_gene_symbol']
					
			gene_id = GENE_SYMBOLS[gene_symbol]
			url = 'http://www.ncbi.nlm.nih.gov/gene/%s' % gene_id

		elif meta['id'].startswith('dz:'):
			do_id = meta.get('do_id', None)
			if do_id is not None:
				url = 'http://disease-ontology.org/term/%s' % do_id
			else:
				url = 'https://www.google.com/search?q=%s' % self.name.replace(' ', '+')
		else:
			db_id = meta.get('drugbank_id', None)
			pubchem_cid = meta.get('pubchem_cid', None)
			if db_id is not None:
				url = 'http://www.drugbank.ca/drugs/%s' % db_id
			elif pubchem_cid is not None:
				url = 'https://pubchem.ncbi.nlm.nih.gov/compound/%s' % pubchem_cid
			else:
				url = 'https://www.google.com/search?q=%s' % self.name.replace(' ', '+')
		return url	



def get_matrix(uids, genes, na_val=0):
	## retrieve a matrix based on uids of signatures and genes
	mat = np.zeros((len(genes), len(uids)))
	projection ={'id':True, '_id':False, 'chdir': True,
		'hs_gene_symbol':True, 'mm_gene_symbol':True, 'organism':True,
		'disease_name':True, 'drug_name':True}

	for j, uid in enumerate(uids):
		sig = DBSignature(uid, projection=projection)
		vals = sig.get_gene_vals(genes, na_val=na_val)
		mat[:, j] = vals
	return mat


def make_download_file(category, format, outfn):
	## to generate files for downloading

	all_sigs = []
	for uid in ALL_UIDS:
		if uid.startswith(category):
			sig = DBSignature(uid) # Signature instance
			if sig.has_chdir():
				sig.fill_top_genes(600)
				all_sigs.append(sig.to_dict(format=format))

	if format == 'gmt':
		with open (outfn, 'w') as out:
			for dict_data in all_sigs:
				line_up = [ dict_data['name'] + '-up', dict_data['id'] ] + map(lambda x:x[0], dict_data['up_genes'])
				out.write('\t'.join(line_up) + '\n')
				line_dn = [ dict_data['name'] + '-dn', dict_data['id'] ] + map(lambda x:x[0], dict_data['down_genes'])
				out.write('\t'.join(line_dn) + '\n')
	else:
		json.dump(all_sigs, open(outfn, 'wb'))
	return len(all_sigs)

def make_all_download_files():
	file_meta = {}
	d_fn = {
		'gene': 'Single_gene_perturbations',
		'dz': 'Disease_signatures',
		'drug': 'Single_drug_perturbations',
	}

	for category in ['gene', 'dz', 'drug']:
		for format in ['json', 'gmt']:
			outfn = 'downloads/%s.%s' % (d_fn[category], format)
			if not os.path.isfile(outfn):
				num_sigs = make_download_file(category, format, outfn)
				print num_sigs
			print category, format, 'finished'
	return


## test
# '''
# gs = DBSignature('gene:24')
# # gs.fill_top_genes(600)
# # print gs.to_json()

# t0 = time.time()

# # d_uid_scores = gs.calc_all_scores(4, 600)
# # print d_uid_scores.items()[:5]
# d_uid_sigs = load_all_db_sigs(nprocs=1)

# tt = time.time()
# # print(len(d_uid_scores))
# print len(d_uid_sigs)
# print('time passed:', tt-t0)

# res = gs.calc_all_scores(d_uid_sigs)
# print res.items()[:4]

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
# print gs.post_to_paea()
# print gs.meta
# print gs.post_to_cds2()

# make_download_file('dz', 'json')
# make_download_file('dz', 'gmt')
# '''



