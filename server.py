## python API for the mongodb
import sys, json
# import numpy as np
# from collections import OrderedDict
from flask import Flask, request
from pymongo import MongoClient

from crossdomain import crossdomain

client = MongoClient('mongodb://127.0.0.1:27017/')
db = client['microtask_signatures']
coll = db['signatures']

app = Flask(__name__)
app.debug = True

@app.route('/microtask_geo_webapp/api', methods=['POST', 'GET'])
@crossdomain(origin='*')
def post_signature():
	if request.method == 'POST':
		data = json.loads(request.data)
		return

	elif request.method == 'GET':
		uid = request.args.get('id', '')
		cutoff = request.args.get('topn', '') # number of genes to return
		if cutoff == '':
			cutoff = 600
		else:
			cutoff = int(cutoff)

		doc = coll.find_one({'id':uid}, {'_id':False})
		if doc is not None:
			chdir = doc['chdir']
			del doc['chdir']
			chdir = sorted(chdir.items(), key=lambda x: abs(x[1]), reverse=True)
			doc['up_genes'] = []
			doc['down_genes'] = []
			for gene, val in chdir[0:cutoff]:
				if val > 0: doc['up_genes'].append( (gene, val) )
				else: doc['down_genes'].append( (gene, val) )

			return json.dumps(doc)
		else: # bad request
			return ('', 400, '')

if __name__ == '__main__':
	if len(sys.argv) > 1:
		port = int(sys.argv[1])
	else:
		port = 5050
	if len(sys.argv) > 2:
		host = sys.argv[2]
	else:
		host = '127.0.0.1'
	app.run(host=host, port=port)
