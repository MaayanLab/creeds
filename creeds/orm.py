'''
ORMs for signature, signatures in the MongoDB and collection of signatures.
'''

import os, sys, json
import time
from itertools import chain
from collections import Counter
import cPickle as pickle
import numpy as np
import pandas as pd
import requests

from .gene_converter import *

## connect to mongodb via pymongo.MongoClient imported from the module
from creeds import conn
COLL = conn['microtask_signatures'].signatures

ALL_UIDS = COLL.find(
	{'$and': [
		{'chdir_sva_exp2': {'$exists': True}}, 
		{'version': {'$in':['1.0', '1.2']}},
		{"incorrect": {"$ne": True}}
	]},
	{'id': True}).distinct('id')

print len(ALL_UIDS)


## load gene symbol to gene ID conversion dict
GENE_SYMBOLS = load_gene_symbol_dict()

## Fields in the mongodb for interal use only
FIELDS_EXCLUDE = ['_id', 
	'limma', 'limma_sva', 'limma_norm', 'limma_combat',
	'fold_changes', 'log2FC_norm',
	'chdir', 'chdir_combat_exp2',
	'pvca', 'pvca_sva', 'pvca_combat']

PROJECTION_EXCLUDE = dict(zip(FIELDS_EXCLUDE, [False] * len(FIELDS_EXCLUDE)))


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
	## assumes l1, l2 are sets of upper case genes
	up = len(l1 & l2)
	dn = len(l1 | l2)
	if dn == 0: # to prevent 0 division error
		return 0
	else:
		return float(up)/dn

def signed_jaccard(s1, s2):
	## signed jaccard index for signatures
	j1 = jaccard(s1._up_genes, s2._up_genes)
	j2 = jaccard(s1._dn_genes, s2._dn_genes)
	j3 = jaccard(s1._dn_genes, s2._up_genes)
	j4 = jaccard(s1._up_genes, s2._dn_genes)
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

def load_and_fill_sigs(uids):
	## retrieve signatures and chdir in batch
	docs = COLL.find({'id': {'$in': uids}}, PROJECTION_EXCLUDE)
	d_uid_sigs = {}
	for doc in docs:
		sig = DBSignature(None, doc=doc)
		sig.fill_top_genes()
		uid = doc['id']
		d_uid_sigs[uid] = sig
	return d_uid_sigs


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
		self._up_genes = set(map(lambda x: x.upper(), up_genes))
		self._dn_genes = set(map(lambda x: x.upper(), dn_genes))

	def calc_all_scores(self, db_sig_collection_iteritems):
		'''
		Calcuated signed jaccard score for this signatures against 
		all DBsignatures in `db_sig_collection_iteritems`, generator from 
		a DBSignatureCollection instance or a list of DBSignatureCollections.
		'''
		uid_scores = []
		for uid, sig in db_sig_collection_iteritems:
			score = signed_jaccard(self, sig)
			uid_scores.append((uid, score))
		return dict(uid_scores)

	def get_query_results(self, db_sig_collection, direction='similar'):
		'''
		Handle querying signatures from the DB with custom up/down genes,
		return a list of objects
		'''
		if type(db_sig_collection) == list:
			db_sig_collection_iteritems = chain(*[dbc.iteritems() for dbc in db_sig_collection])
		else:
			db_sig_collection_iteritems = db_sig_collection.iteritems()
		

		d_uid_score = self.calc_all_scores(db_sig_collection_iteritems)

		scores = np.array(d_uid_score.values())
		uids = np.array(d_uid_score.keys())
		uid_data = [] # a list of meta data {} sorted by score

		# mask for signs of scores
		if direction == 'similar':
			score_sign_mask = scores > 0
		elif direction == 'opposite':
			score_sign_mask = scores < 0
		# sort uids by abs(scores) in descending order
		srt_idx = np.abs(scores[score_sign_mask]).argsort()[::-1]
		scores = scores[score_sign_mask][srt_idx]
		uids = uids[score_sign_mask][srt_idx]

		# retrieve meta-data for all uids
		projection ={'geo_id':True, 'id':True, '_id':False,
			'hs_gene_symbol':True, 'mm_gene_symbol':True, 'organism':True, 
			'disease_name':True, 'drug_name':True, 'do_id':True,
			'drugbank_id':True, 'pubchem_cid':True}
		uid_docs = COLL.find({'id': {'$in': uids.tolist()}}, projection)

		uid_docs = list(uid_docs)
		# make uid_docs have the same order of id with uids
		uids = uids.tolist()
		uid_docs_ = [None] * len(uid_docs)
		for uid_doc in uid_docs:
			idx = uids.index(uid_doc['id'])
			uid_docs_[idx] = uid_doc
		uid_docs = uid_docs_

		for doc, score in zip(uid_docs, scores):
			sig_ = DBSignature(None, doc=doc)
			meta = {
				'id': sig_.meta['id'],
				'geo_id': sig_.meta['geo_id'],
				'name': [sig_.name, sig_.get_url()], # [name, url]
				'signed_jaccard': float('%.5f'%score)
			}
			uid_data.append(meta)
		return uid_data


class DBSignature(Signature):
	'''
	Signature instance from the mongodb.
	'''
	chdir_field = 'chdir_sva_exp2'

	def __init__(self, uid, projection=PROJECTION_EXCLUDE, doc=None):
		## the constructor also act as a way to query mongodb using
		## the id and return desirable fields by specifying projection
		if doc is None: ## if doc is given, do not retrieve from DB
			doc = COLL.find_one({'id':uid}, projection)
		name = find_name(doc)
		if self.chdir_field in doc:
			chdir = doc[self.chdir_field]
			del doc[self.chdir_field]
			self.chdir = chdir
		Signature.__init__(self, name, doc)

	def has_chdir(self):
		## assert if a signature has the chdir field
		if hasattr(self, 'chdir'): return True
		else: return False

	def is_filled(self):
		## assert if top genes are filled
		return len(self.up_genes) > 0 and len(self.dn_genes) > 0

	def fill_top_genes(self, cutoff=600):
		## get top up/dn genes based on a rank cutoff
		if not self.is_filled():
			for gene, val in zip(self.chdir['genes'], self.chdir['vals'])[:cutoff]:
				if val > 0: 
					self.up_genes.append( (gene, val) )
					self._up_genes.add( gene.upper() )
				else: 
					self.dn_genes.append( (gene, val) )
					self._dn_genes.add( gene.upper() )


	def to_json(self, meta_only=False):
		## to export the document into json
		json_data = self.meta
		if not meta_only:
			self.fill_top_genes()
			json_data['up_genes'] = self.up_genes
			json_data['down_genes'] = self.dn_genes
		return json.dumps(json_data)

	def to_dict(self, format='gmt'):
		## method to generate files for downloading
		if format == 'gmt':
			dict_data = {'name': self.name, 'id': self.meta['id']}			
		else:
			dict_data = self.meta

		self.fill_top_genes()
		dict_data['up_genes'] = self.up_genes
		dict_data['down_genes'] = self.dn_genes
		return dict_data


	def calc_all_scores(self, db_sig_collection, cutoff=600):
		## calcuated signed jaccard score for this signatures against 
		## all signatures in the database
		self.fill_top_genes(cutoff)
		uid_scores = []
		for uid, sig in db_sig_collection.items():
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
					
			gene_id = GENE_SYMBOLS.get(gene_symbol, '')
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


class DBSignatureCollection(dict):
	'''
	A collection of DBSignature from the mongodb
	'''
	formats = ['csv', 'json', 'gmt']
	category2name = {
		'gene': 'Single_gene_perturbations',
		'dz': 'Disease_signatures',
		'drug': 'Single_drug_perturbations',
		}

	outfn_path = os.path.dirname(os.path.realpath(__file__)) + '/static/downloads/'

	def __init__(self, filter_, name):
		'''
		`filter_` should be a mongo query
		'''
		self.filter_ = filter_
		self.name = name

		cur = COLL.find(self.filter_, PROJECTION_EXCLUDE)
		uids = cur.distinct('id')
		# Load signatures 
		for doc in cur:
			sig = DBSignature(None, doc=doc)
			sig.fill_top_genes()
			uid = doc['id']
			self[uid] = sig

		categories = map(lambda x:x.split(':')[0], self.keys())
		self.categories = set(categories)
		self.category_count = dict(Counter(categories))
		self.get_download_file_meta()

	def get_download_file_meta(self):
		'''to store metadata for generated download files
		'''
		self.download_file_meta = []
		for category, num_sigs in self.category_count.items():
			category_name = self.category2name[category]

			filenames = {format: '%s-%s.%s' % (category_name, 
				self.name, format) for format in self.formats}
			meta = {
				'name': self.name,
				'category': category_name,
				'filenames': filenames,
				'num_sigs': num_sigs
				}

			self.download_file_meta.append(meta)


	def make_download_file(self, category, format, outfn):
		'''to generate files for downloading
		'''
		sigs_this_category = [sig.to_dict(format=format) for uid, sig in self.items() if uid.startswith(category)]
		if format == 'gmt':
			with open (outfn, 'w') as out:
				for dict_data in sigs_this_category:
					line_up = [ dict_data['name'] + '-up', dict_data['id'] ] + map(lambda x:x[0], dict_data['up_genes'])
					out.write('\t'.join(line_up) + '\n')
					line_dn = [ dict_data['name'] + '-dn', dict_data['id'] ] + map(lambda x:x[0], dict_data['down_genes'])
					out.write('\t'.join(line_dn) + '\n')
		
		elif format == 'csv': # annotations only
			df = pd.DataFrame.from_records(sigs_this_category)\
				.drop(['up_genes', 'down_genes'], axis=1)\
				.set_index('id')
			df['pert_ids'] = df['pert_ids'].map(lambda x: ','.join(x))
			df['ctrl_ids'] = df['ctrl_ids'].map(lambda x: ','.join(x))
			df.to_csv(outfn, encoding='utf-8')

		else:
			json.dump(sigs_this_category, open(outfn, 'wb'))
		return (outfn, len(sigs_this_category))

	def make_all_download_files(self):
		'''to generate all 3 formats of files for each category
		'''
		for category in self.categories:
			for format in self.formats:
				outfn = '%s/%s-%s.%s' % (self.outfn_path, 
					self.category2name[category], self.name, format)
				if not os.path.isfile(outfn):
					(filename, num_sigs) = self.make_download_file(category, format, outfn)
					print filename, num_sigs
				print outfn, 'finished'
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
# d_cat_names = make_autocomplete()
# json.dump(d_cat_names, open('data/autoCompleteList.json', 'wb'))
# '''
