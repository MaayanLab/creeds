'''
ORMs for signature, signatures in the MongoDB and collection of signatures.
'''
import os, sys, json
from collections import Counter

import numpy as np
import pandas as pd
import scipy.sparse as sp
import requests
from joblib import Parallel, delayed

from .gene_converter import *
from .matrix_ops import (fast_jaccard, fast_signed_jaccard)

## connect to mongodb via pymongo.MongoClient imported from the module
from creeds import conn

################################ Global variables ################################
COLL = conn['microtask_signatures'].signatures
COLL_GENES = conn['microtask_signatures'].genes

ALL_GENES = COLL_GENES.find_one({'case_sensitive': {'$exists':True}})['case_sensitive']
ALL_GENES = np.array(ALL_GENES)
ALL_GENES_I = COLL_GENES.find_one({'case_insensitive': {'$exists':True}})['case_insensitive']
ALL_GENES_I = np.array(ALL_GENES_I)

ALL_UIDS = COLL.find(
	{'$and': [
		{'chdir_sva_exp2': {'$exists': True}}, 
		{'version': {'$in':['1.0', '1.1', '1.2', '2.0']}},
		{"incorrect": {"$ne": True}}
	]},
	{'id': True}).distinct('id')

print '# signatures: %d, # genes(s): %d, # genes(i) : %d' \
	% (len(ALL_UIDS), len(ALL_GENES), len(ALL_GENES_I))


## load gene symbol to gene ID conversion dict
GENE_SYMBOLS = load_gene_symbol_dict()

## Fields in the mongodb for interal use only
FIELDS_EXCLUDE = ['_id', 
	'limma', 'limma_sva', 'limma_norm', 'limma_combat',
	'fold_changes', 'log2FC_norm',
	'chdir', 'chdir_combat_exp2',
	'pvca', 'pvca_sva', 'pvca_combat']

PROJECTION_EXCLUDE = dict(zip(FIELDS_EXCLUDE, [False] * len(FIELDS_EXCLUDE)))

################################ Util functions ################################
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
			name = doc.get('mm_gene_symbol', None)
			if name is None:
				name = doc['hs_gene_symbol']
	elif prefix == 'dz':
		name = doc['disease_name']
	else:
		name = doc['drug_name']

	if type(name) == list: # predicted signatures
		# extract name fields and convert to string
		name = [item['name'] for item in name]
	return name


################################ Classes ################################
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


	def init_vectors(self):	
		'''Init binary vectors representing of the siganture, 
		for fast computation of jaccard.'''
		v_up = np.zeros(len(ALL_GENES_I), dtype=np.int8)
		up_genes_i = map(lambda x: x.upper(), self.up_genes)
		v_up[np.in1d(ALL_GENES_I, up_genes_i)] = 1

		v_dn = np.zeros(len(ALL_GENES_I), dtype=np.int8)
		dn_genes_i = map(lambda x: x.upper(), self.dn_genes)
		v_dn[np.in1d(ALL_GENES_I, dn_genes_i)] = 1

		self.v_up = v_up
		self.v_dn = v_dn

	def calc_all_scores(self, db_sig_collection):
		'''
		Calcuated signed jaccard score for this signatures against 
		a DBSignatureCollection instance or a list of DBSignatureCollection instances.
		'''
		uid_scores = []
		if type(db_sig_collection) != list: 
			# a single DBSignatureCollection instance
			scores = fast_signed_jaccard(db_sig_collection.mat_up, db_sig_collection.mat_dn, 
				self.v_up, self.v_dn)
			uids = db_sig_collection.uids

		else:
			# a list of DBSignatureCollection instances
			# stack sparse matrices first
			mat_up = sp.vstack([dbsc.mat_up for dbsc in db_sig_collection])
			mat_dn = sp.vstack([dbsc.mat_dn for dbsc in db_sig_collection])
			scores = fast_signed_jaccard(mat_up, mat_dn, 
				self.v_up, self.v_dn)
			uids = []
			for dbsc in db_sig_collection:
				uids.extend(dbsc.uids)
		
		uid_scores = zip(uids, scores)
		return dict(uid_scores)

	def get_query_results(self, db_sig_collection, direction='similar'):
		'''
		Handle querying signatures from the DB with custom up/down genes,
		return a list of objects
		'''

		d_uid_score = self.calc_all_scores(db_sig_collection)

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

	def get_vector_indexes(self, cutoff=600):
		'''
		Get indexes of up and down genes in ALL_GENES_I
		Used for construction of sparse matices in DBSignatureCollection.mat_*
		'''	
		genes_i = np.array(map(lambda x:x.upper(), self.chdir['genes'][:cutoff]))
		vals = np.array(self.chdir['vals'][:cutoff])
		up_genes_i = genes_i[vals > 0]
		dn_genes_i = genes_i[vals < 0]
		return np.in1d(ALL_GENES_I, up_genes_i), np.in1d(ALL_GENES_I, dn_genes_i)

	def init_cs_vectors(self, cutoff=600):
		'''Init case sensitive vectors with CD values. 
		This vector is intended to for exporting purpose used in to_json, to_dict.
		'''
		if not hasattr(self, 'v_cs'):
			v_cs = np.zeros(len(ALL_GENES), dtype=np.float32)
			genes = self.chdir['genes'][:cutoff]
			genes, uniq_idx = np.unique(genes, return_index=True)
			vals = np.array(self.chdir['vals'][:cutoff])[uniq_idx]

			v_cs[np.in1d(ALL_GENES, genes)] = vals
			self.v_cs = sp.lil_matrix(v_cs)

	def fill_top_genes(self):
		'''Get top up/dn genes from `v_cs`
		'''
		# get mask of non zero index
		mask_non_zero = (self.v_cs != 0).toarray().ravel()
		# retrieve CD vals and genes
		vals = self.v_cs[0, mask_non_zero].toarray().ravel()
		genes = ALL_GENES[mask_non_zero]
		# sort CD vals on abs(vals)
		srt_idx = np.abs(vals).argsort()[::-1]

		up_genes = []
		dn_genes = []
		for gene, val in zip(genes[srt_idx].tolist(), vals[srt_idx].tolist()):
			if val > 0: 
				up_genes.append( (gene, val) )
			else: 
				dn_genes.append( (gene, val) )
		return up_genes, dn_genes


	def clear(self, cutoff=600):
		'''Clear unnecessary fields to reduce RAM usage
		'''
		self.init_cs_vectors(cutoff=cutoff)
		del self.chdir

	def to_json(self, meta_only=False):
		## to export the document into json
		json_data = self.meta
		if not meta_only:
			up_genes, dn_genes = self.fill_top_genes()
			json_data['up_genes'] = up_genes
			json_data['down_genes'] = dn_genes

		return json.dumps(json_data)

	def to_dict(self, format='gmt'):
		## method to generate files for downloading
		if format == 'gmt':
			dict_data = {'name': self.name, 'id': self.meta['id']}			
		else:
			dict_data = self.meta

		up_genes, dn_genes = self.fill_top_genes()
		dict_data['up_genes'] = up_genes
		dict_data['down_genes'] = dn_genes
		return dict_data


	def post_to_paea(self, cutoff=2000):
		## post top n genes to PAEA and return a PAEA url
		## return None if instance has no chdir
		post_url = 'http://amp.pharm.mssm.edu/Enrichr/addList'
		base_url = 'http://amp.pharm.mssm.edu/PAEA?id='
		paea_url = None
		if self.has_chdir():
			up_genes, dn_genes = self.fill_top_genes()
			gene_list = []
			for gene, coef in up_genes + dn_genes:
				gene_list.append( '%s,%s\n'% (gene, coef) )
			gene_list = ''.join(gene_list)
			data = {'list': gene_list, 'inputMethod': "PAEA", 'description': self.name}
			r = requests.post(post_url, files=data)
			paea_url = base_url + str(json.loads(r.text)['userListId'])
		return paea_url

	def post_to_cds2(self, cutoff=2000):
		## post top n genes to L1000CDS2 API and return a CDS2 url
		url = 'http://amp.pharm.mssm.edu/L1000CDS2/query'
		cds2_url = None
		if self.has_chdir():
			up_genes, dn_genes = self.fill_top_genes()
			data = {
				"genes": map(lambda x: x[0].upper(), up_genes + dn_genes), 
				"vals":  map(lambda x: x[1], up_genes + dn_genes)
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
		uid = meta['id']
		if ':P' not in uid: # not v2.0 signature
			if uid.startswith('gene:'):
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

			elif uid.startswith('dz:'):
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

		else: # v2.0 signature
			if uid.startswith('gene:'): key = 'hs_gene_symbol'
			elif uid.startswith('dz:'): key = 'disease_name'
			else: key = 'drug_name'
			# becas concept ids
			cids = [':'.join(item['cid'].split(':')[:2]) for item in meta[key]]
			url_root = 'http://bioinformatics.ua.pt/becas/api/concept/redirect/'
			url = [url_root + cid for cid in cids]

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

	def __init__(self, filter_, name, limit=None):
		'''
		`filter_` should be a mongo query
		'''
		self.filter_ = filter_
		self.name = name

		if not limit:
			cur = COLL.find(self.filter_, PROJECTION_EXCLUDE)
		else:
			cur = COLL.find(self.filter_, PROJECTION_EXCLUDE).limit(limit)
		# to preserve orders
		self.uids = cur.distinct('id')
		
		# sparse matrices
		sparse_mat_shape = (len(self.uids), len(ALL_GENES_I))
		mat_up = sp.lil_matrix(sparse_mat_shape, dtype=np.int8)
		mat_dn = sp.lil_matrix(sparse_mat_shape, dtype=np.int8)

		# Load signatures 
		for i, doc in enumerate(cur):
			sig = DBSignature(None, doc=doc)
			# sig.fill_top_genes()
			# fill the sparse matrices
			up_idx, dn_idx = sig.get_vector_indexes()
			mat_up[i, up_idx] = 1
			mat_dn[i, dn_idx] = 1
			# clear `chdir` field and add `v_cs` for exporting
			sig.clear(cutoff=600)

			uid = doc['id']
			self[uid] = sig

		# convert to CSR format for fast compuatation
		self.mat_up = mat_up.tocsr()
		self.mat_dn = mat_dn.tocsr()
		del mat_up, mat_dn

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
		return


	def make_download_file(self, category, format, outfn):
		'''to generate files for downloading
		'''
		sigs_this_category = [sig.to_dict(format=format) for uid, sig in self.items() if uid.startswith(category)]
		if format == 'gmt':
			with open (outfn, 'w') as out:
				for dict_data in sigs_this_category:
					if type(dict_data['name']) == list: # join list to string for v2.0
						term_up = '|'.join(dict_data['name']) + '-up'
						term_dn = '|'.join(dict_data['name']) + '-dn'
					else:
						term_up = dict_data['name'] + '-up'
						term_dn = dict_data['name'] + '-dn'

					line_up = [ term_up, dict_data['id'] ] + map(lambda x:x[0], dict_data['up_genes'])
					line_dn = [ term_dn, dict_data['id'] ] + map(lambda x:x[0], dict_data['down_genes'])
					out.write('\t'.join(line_up) + '\n')
					out.write('\t'.join(line_dn) + '\n')
		
		elif format == 'csv': # annotations only
			if 'entities' in sigs_this_category[0]:
				# predicted signatures
				# need to unpack fields that are dict
				sigs_this_category_ = []
				for doc in sigs_this_category:
					del doc['entities']
					for field in ['hs_gene_symbol', 'disease_name', 'drug_name', 'cell_type']:
						if field in doc:
							doc[field] = '|'.join([item['name'] for item in doc[field]])
					sigs_this_category_.append(doc)
				sigs_this_category = sigs_this_category_ 

			df = pd.DataFrame.from_records(sigs_this_category)\
				.drop(['up_genes', 'down_genes'], axis=1)\
				.set_index('id')
			df['pert_ids'] = df['pert_ids'].map(lambda x: '|'.join(x))
			df['ctrl_ids'] = df['ctrl_ids'].map(lambda x: '|'.join(x))
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


