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
