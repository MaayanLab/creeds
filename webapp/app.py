## python API for the mongodb
import os, sys, json
# import numpy as np
# from collections import OrderedDict
from flask import Flask, request

from orm_utils import *
from crossdomain import crossdomain


app = Flask(__name__, static_url_path='', static_folder=os.getcwd())
app.debug = True

@app.route('/')
def root():
	return app.send_static_file('index.html')


@app.route('/api', methods=['GET'])
@crossdomain(origin='*')
def retrieve_signature():
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
			return sig.to_json(meta_only=False)
		else: # bad request
			return ('', 400, '')


@app.route('/autoCompleteList', methods=['GET', 'POST'])
@crossdomain(origin='*')
def get_all_names():
	if request.method == 'GET':
		## get a list of names for signatures in mongodb with chdir
		names = []
		for uid in ALL_UIDS:
			sig = DBSignature(uid)
			if sig.has_chdir():
				names.append(sig.name)
		return json.dumps(names)


@app.route('/search', methods=['GET', 'POST'])
@crossdomain(origin='*')
def search():
	## handle searching signatures with either custom genes or a document in the db
	def get_search_results(sig):
		d_uid_score = sig.calc_all_scores()
		uid_data = [] # a list of meta data {} sorted by score
		for uid, score in sorted(d_uid_score.items(), key=lambda x:x[1]): # small to large signed_jaccard
			projection ={'geo_id':True, 'id':True, '_id':False,
				'hs_gene_symbol':True, 'disease_name':True, 'drug_name':True}
			sig_ = DBSignature(uid, projection=projection)
			meta = sig_.meta
			meta['name'] = sig_.name
			meta['signed_jaccard'] = score
			uid_data.append(meta)
		return uid_data

	if request.method == 'GET': # search using an id in the mongodb
		uid = request.args.get('id', '')
		if uid in ALL_UIDS:
			sig = DBSignature(uid) # Signature instance
		else: # invalid uid
			sig = None

	elif request.method == 'POST': # search using custom up/dn gene list
		data = json.loads(request.data)
		up_genes = data['up_genes']
		dn_genes = data['dn_genes']
		name = data['name']
		meta = data['meta']
		sig = Signature(name, meta, up_genes, dn_genes)

	if sig is not None:
		uid_data = get_search_results(sig)
		return json.dumps(uid_data)
	else:
		return ('', 400, '')



if __name__ == '__main__':
	if len(sys.argv) > 1:
		port = int(sys.argv[1])
	else:
		port = 5000
	if len(sys.argv) > 2:
		host = sys.argv[2]
	else:
		host = '127.0.0.1'
		# host = '0.0.0.0'
	app.run(host=host, port=port)
