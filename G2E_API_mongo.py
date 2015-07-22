## hit the G2E API with the meta data from mongodb microtask_signatures

import os, sys
import numpy as np
import time
import urllib2
import json
from pymongo import MongoClient
import cPickle as pickle
from pprint import pprint


DATADIR = 'output/microtask_dz_jsons_limma_norm' 
PREFIX = 'dz:'
BASE_URL = 'http://127.0.0.1:8086/g2e/full?'

client = MongoClient('mongodb://127.0.0.1:27017/')
db = client['microtask_signatures']
coll = db['signatures']
all_uid = coll.distinct('id')


def get_json(doc):
	url = BASE_URL
	## required args
	url += 'accession=%s&'%doc['geo_id']
	url += 'platform=%s&'%doc['platform']
	url += 'control=%s&'%'-'.join(doc['ctrl_ids'])
	url += 'experimental=%s'%'-'.join(doc['pert_ids'])
	url += '&cutoff=None' ## new args to retrieve the full chdir
	print url
	response = urllib2.urlopen(url)
	data = json.load(response)
	# print 'API status:', data['status']
	return data


uids = []
projection = {'_id':False, 'chdir':True, 'geo_id':True, 'ctrl_ids':True, 'pert_ids':True,'platform':True}
for uid in all_uid:

	if uid.startswith(PREFIX):
		uids.append( uid)
		id = uid.split(':')[1]
		json_fn = DATADIR + '/' + id+'.json'
		doc = coll.find_one({'id': uid}, projection)
		if 'chdir' in doc:
			try:
				json_data = get_json(doc)
				json.dump(json_data, open(json_fn,'wb'))
				print json_fn, 'success'
			except Exception, e:
				print uid, 'error', e
				pass
			



# print len(uids)
# print uids[-5:]

