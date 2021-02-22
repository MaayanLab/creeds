## python API for the mongodb
import os, sys, json, pickle
# Reset timezone to EDT 
import time
os.environ['TZ'] = 'US/Eastern'
time.tzset()

import logging
from collections import OrderedDict

import clustergram
from .crossdomain import crossdomain

import re
from werkzeug.routing import Rule, RequestRedirect

## super hacky way to make the URL of the app case insensitive...
## http://librelist.com/browser/flask/2011/6/24/case-insensitive-routing/#198dd20c7198760b3e2f5d5ada19b7f9
class CIRule(Rule):
	def compile(self):
		Rule.compile(self)
		self._regex = re.compile(self._regex.pattern, 
			re.UNICODE | re.IGNORECASE)

from flask import (Flask, request, Response, send_from_directory)

class CIFlask(Flask):
    url_rule_class = CIRule


SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))

from pymongo import MongoClient

ENTER_POINT = '/CREEDS'
app = CIFlask(__name__, static_url_path=ENTER_POINT, static_folder='static')

# Get config from object
app.config.from_object(os.environ['CONFIG_OBJ'])
# Make connection with MongoDB
conn = MongoClient(app.config['DATABASE_URI'])
# Import models and utils
from .orm import *
from .utils import *

app.logger.setLevel(logging.INFO)

if not app.debug:
	# Add a logging handler for production to log out INFOs
	logging_handler = logging.StreamHandler()
	formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
	logging_handler.setFormatter(formatter)
	app.logger.addHandler(logging_handler)
	
app.logger.info('Loading DBSignatureCollections')
with open(SCRIPT_PATH+'/static/downloads/d_dbsc.pickle', 'r') as fr:
	d_dbsc = pickle.load(fr)
app.logger.info('DBSignatureCollections loaded')

def create_globals():
	''' This function can be used to construct the data
	'''
	# Load globals DBSignatureCollection instances
	global d_dbsc
	d_dbsc = OrderedDict() # {name : DBSignatureCollection instance}
	
	app.logger.info('# signatures: %d, # genes(s): %d, # genes(i) : %d' \
		% (len(ALL_UIDS), len(ALL_GENES), len(ALL_GENES_I)))

	app.logger.info('Loading DBSignatureCollections')
	app.logger.info('\t'.join(['Collection name','Number of signatures','Size(MB)']))

	for params in app.config['DBSC_PARAMS']:
		collection_name = params['name']
		d_dbsc[collection_name] = DBSignatureCollection(**params)
		app.logger.info('%s\t%d\t%.2f' % (collection_name, 
			len(d_dbsc[collection_name]), 
			sys.getsizeof(d_dbsc[collection_name])/1e6))

	app.logger.info('DBSignatureCollections loaded')
	
	if app.config['MAKE_DOWNLOAD_FILES']:
		for dbsc in d_dbsc.values():
			dbsc.make_all_download_files()

	with open(SCRIPT_PATH+'/static/downloads/d_dbsc.pickle', 'w') as fw:
		pickle.dump(d_dbsc, fw)

	return


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
			sig.init_cs_vectors(cutoff)
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
	## handle searching signatures using string and query signatures with up/down genes 
	if request.method == 'GET': # search signatures using string
		search_string = request.args.get('q', '')
		search_dict = {'$and': [
			{'chdir_sva_exp2': {'$exists': True}},
			{"incorrect": {"$ne": True}},
			{"$or":[
				{'hs_gene_symbol' : {"$regex": search_string, "$options":"i"}},
				{'hs_gene_symbol' : {'$elemMatch': {'name': {"$regex": search_string, "$options":"i"}}}},
				{'mm_gene_symbol' : {"$regex": search_string, "$options":"i"}},
				{'disease_name' : {"$regex": search_string, "$options":"i"}},
				{'disease_name' : {'$elemMatch': {'name': {"$regex": search_string, "$options":"i"}}}},
				{'drug_name' : {"$regex": search_string, "$options":"i"}},
				{'drug_name' : {'$elemMatch': {'name': {"$regex": search_string, "$options":"i"}}}},
				]}
		]}
		
		docs = []
		for doc in COLL.find(search_dict, {'id':True,'_id':False}):
			uid = doc['id']
			sig = DBSignature(uid) # Signature instance
			docs.append(sig.meta)

		return Response(json.dumps(docs), mimetype='application/json')

	elif request.method == 'POST': # search using custom up/dn gene list
		data = json.loads(request.data)
		up_genes = data['up_genes']
		dn_genes = data['dn_genes']

		name = data.get('name', None)
		meta = data.get('meta', None)
		direction = data.get('direction', 'similar')
		db_version = data.get('db_version', 'v1.0')
		query_params = dict(direction=direction, db_version=db_version)
		client = data.get('client', 'api')

		sig = Signature(name, meta, up_genes, dn_genes, query_params)
		if sig is None:
			return ('', 400, '')
		else:
			# Save the user signature
			h = sig.save()

			if client == 'api': # perform the query only when client is api
				sig.init_vectors()
				uid_data = sig.get_query_results(d_dbsc)

				return Response(json.dumps(uid_data), mimetype='application/json')

			else: # return the hash of the signature to front-end
				return Response(json.dumps(h), mimetype='application/json')


@app.route(ENTER_POINT + '/result', methods=['GET'])
@crossdomain(origin='*')
def result():
	## retrieve a Signature using a hash, perform query, then return the result
	if request.method == 'GET':
		h = request.args.get('id', None)
		sig = Signature.from_hash(h)
		if sig is None:
			return ('', 400, '')
		else:
			sig.init_vectors()
			uid_data = sig.get_query_results(d_dbsc)
			return Response(json.dumps(uid_data), mimetype='application/json')


@app.route(ENTER_POINT + '/download', methods=['GET'])
def send_download_meta():
	## send meta data for files to download
	all_download_file_meta = []
	for dbsc in d_dbsc.values():
		all_download_file_meta.extend(dbsc.download_file_meta)
	return json.dumps(all_download_file_meta)


@app.route(ENTER_POINT + '/download/<path:filename>')
def download_file(filename):
	## safely send the file for download
	return send_from_directory(SCRIPT_PATH+'/static/downloads', 
		filename, as_attachment=True)


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
			sig.init_cs_vectors(cutoff=2000)
			if app_name == 'paea':
				url = sig.post_to_paea(cutoff=2000)
			elif app_name == 'cds2':
				url = sig.post_to_cds2(cutoff=2000)
		return json.dumps(url)	


@app.errorhandler(500)
def internal_error(exception):
	app.logger.error(exception)

