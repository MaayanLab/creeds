## python API for the mongodb
import os, sys, json

from orm_utils import *
import clustergram
from crossdomain import crossdomain

import re
from werkzeug.routing import Rule, RequestRedirect

## super hacky way to make the URL of the app case insensitive...
## http://librelist.com/browser/flask/2011/6/24/case-insensitive-routing/#198dd20c7198760b3e2f5d5ada19b7f9
class CIRule(Rule):
	def compile(self):
		Rule.compile(self)
		self._regex = re.compile(self._regex.pattern, 
			re.UNICODE | re.IGNORECASE)

from flask import (Flask, request, Response)

class CIFlask(Flask):
    url_rule_class = CIRule


ENTER_POINT = '/CREEDS'
app = CIFlask(__name__, static_url_path=ENTER_POINT, static_folder='static')
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 6


@app.route(ENTER_POINT + '/')
def root():
	return app.send_static_file('index.html')

@app.route(ENTER_POINT + '/api', methods=['GET'])
@crossdomain(origin='*')
def retrieve_signature():
	## to retrieve data and meta data given id of signature like `gene:24`
	if request.method == 'GET':
		uid = request.args.get('id', '')
		cutoff = request.args.get('topn', '') # number of genes to return
		if cutoff == '':
			cutoff = 600
		else:
			cutoff = int(cutoff)

		if uid in ALL_UIDS:
			sig = DBSignature(uid) # Signature instance
			sig.fill_top_genes(cutoff)
			return Response(sig.to_json(meta_only=False), mimetype='application/json')
		else: # bad request
			return ('', 400, '')


@app.route(ENTER_POINT + '/autoCompleteList', methods=['GET', 'POST'])
@crossdomain(origin='*')
def get_all_names():
	## get a list of names for signatures in mongodb with chdir
	if request.method == 'GET':
		d_cat_names = make_autocomplete()
		return json.dumps(d_cat_names)


@app.route(ENTER_POINT + '/search', methods=['GET', 'POST'])
@crossdomain(origin='*')
def search():
	## handle searching signatures with either custom genes or a document in the db
	def get_search_results(sig, direction='similar'):
		d_uid_score = sig.calc_all_scores(d_uid_sigs)
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

	if request.method == 'GET': # search signatures using string
		search_string = request.args.get('q', '')
		search_dict = {
			"$or":[
				{'hs_gene_symbol' : {"$regex": search_string, "$options":"i"}},
				{'mm_gene_symbol' : {"$regex": search_string, "$options":"i"}},
				{'disease_name' : {"$regex": search_string, "$options":"i"}},
				{'drug_name' : {"$regex": search_string, "$options":"i"}}
				]
			}
		docs = []
		for doc in COLL.find(search_dict, {'id':True,'_id':False}):
			uid = doc['id']
			sig = DBSignature(uid) # Signature instance
			if sig.has_chdir():
				docs.append(sig.meta)

		return Response(json.dumps(docs), mimetype='application/json')

	elif request.method == 'POST': # search using custom up/dn gene list
		data = json.loads(request.data)
		up_genes = map(lambda x : x.upper(), data['up_genes'])
		dn_genes = map(lambda x : x.upper(), data['dn_genes'])
		name = data.get('name', None)
		meta = data.get('meta', None)
		direction = data.get('direction', 'similar')

		sig = Signature(name, meta, up_genes, dn_genes)
		uid_data = get_search_results(sig, direction=direction)

		if sig is not None:
			return Response(json.dumps(uid_data), mimetype='application/json')
		else:
			return ('', 400, '')


@app.route(ENTER_POINT + '/geneSigClustergram', methods=['POST'])
@crossdomain(origin='*')
def make_gene_sig_clustergram():
	if request.method == 'POST':
		post_data = json.loads(request.data)
		uids = post_data.get('ids', '')
		genes = post_data.get('genes', '')
		genes = map(lambda x: x.upper(), genes)
		na_val = post_data.get('na_val', 0)
		mat = get_matrix(uids, genes, na_val=na_val)

		json_data = clustergram.clustergram(mat, genes, uids)
		return json.dumps(json_data)


# @app.route('/sigSigClustergram', methods=['GET', 'POST'])
# @crossdomain(origin='*')
## this one shoud probably be pre-computed
@app.route(ENTER_POINT + '/appUrl', methods=['GET'])
@crossdomain(origin='*')
def get_link():
	## to get L1000CDS2 and PAEA url for a DBSignature instance in DB 
	if request.method == 'GET':
		uid = request.args.get('id', '')
		app_name = request.args.get('app', '')
		url = None
		if uid in ALL_UIDS:
			sig = DBSignature(uid) # Signature instance
			if app_name == 'paea':
				url = sig.post_to_paea(cutoff=2000)
			elif app_name == 'cds2':
				url = sig.post_to_cds2(cutoff=2000)
		return json.dumps(url)	



if __name__ == '__main__':
	if len(sys.argv) > 1:
		port = int(sys.argv[1])
	else:
		port = 5000
	if len(sys.argv) > 2:
		host = sys.argv[2]
	else:
		host = '127.0.0.1'

	if len(sys.argv) > 3:
		mode = sys.argv[3]
	else:
		mode = 'dev'

	if mode == 'dev':
		app.debug = True
	else:
		app.debug = False	

	make_all_download_files()
	global d_uid_sigs
	d_uid_sigs = load_all_db_sigs(nprocs=1)
	print 'd_uid_sigs loaded,', len(d_uid_sigs)

	app.run(host=host, port=port)

