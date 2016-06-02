import os, json
import unittest
import pandas as pd
from flask import (request, Response)
# Assumes env var for config is set
from creeds import app

ENTRY_POINT = '/creeds'

UP_GENES = ['KIAA0907','KDM5A','CDC25A','EGR1','GADD45B','RELB','TERF2IP','SMNDC1','TICAM1','NFKB2',
	'RGS2','NCOA3','ICAM1','TEX10','CNOT4','ARID4B','CLPX','CHIC2','CXCL2','FBXO11','MTF2',
	'CDK2','DNTTIP2','GADD45A','GOLT1B','POLR2K','NFKBIE','GABPB1','ECD','PHKG2','RAD9A',
	'NET1','KIAA0753','EZH2','NRAS','ATP6V0B','CDK7','CCNH','SENP6','TIPARP','FOS','ARPP19',
	'TFAP2A','KDM5B','NPC1','TP53BP2','NUSAP1']

DN_GENES = ['SCCPDH','KIF20A','FZD7','USP22','PIP4K2B','CRYZ','GNB5','EIF4EBP1','PHGDH','RRAGA',
	'SLC25A46','RPA1','HADH','DAG1','RPIA','P4HA2','MACF1','TMEM97','MPZL1','PSMG1','PLK1',
	'SLC37A4','GLRX','CBR3','PRSS23','NUDCD3','CDC20','KIAA0528','NIPSNAP1','TRAM2','STUB1',
	'DERA','MTHFD2','BLVRA','IARS2','LIPA','PGM1','CNDP2','BNIP3','CTSL1','CDC25B','HSPA8',
	'EPRS','PAX8','SACM1L','HOXA5','TLE1','PYGL','TUBB6','LOXL1']

class TestApiRequests(unittest.TestCase):
	def test_get_request_args(self):
		with app.test_request_context(ENTRY_POINT +'/api?id=gene:27'):
			self.assertEquals(request.path, ENTRY_POINT + '/api')
			self.assertEquals(request.args['id'], 'gene:27')



class TestStringSearch(unittest.TestCase):
	'''
	Test the API of search the DB using string
	'''
	def setUp(self):
		self.app = app.test_client()

	def test_search(self):
		# Send GET request and get Response object
		resp = self.app.get(ENTRY_POINT + '/search?q=TP53')
		
		self.assertEquals(resp.status_code, 200)
		self.assertEquals(resp.mimetype, 'application/json')

		# Parse the repsonse data
		signatures = json.loads(resp.data.decode())
		signature1 = signatures[0]
		self.assertEquals(type(signatures), list)
		self.assertEquals(type(signature1), dict)

		# Get shared keys in all signatures
		shared_keys = reduce(lambda sig1, sig2: set(sig1) & set(sig2), signatures)
		expected_keys = set(['id', 'geo_id', 'pert_ids', 'ctrl_ids', 'platform', 'version'])
		self.assertTrue(shared_keys.issuperset(expected_keys))


class TestQueryUsingGeneLists(unittest.TestCase):
	'''
	Test the API handling query of DB using up/down gene lists
	'''

	payload = {
		'up_genes': UP_GENES,
		'dn_genes': DN_GENES,
		'direction': 'opposite',
		'db_version': 'v1.0'
		}
	
	def setUp(self):
		self.app = app.test_client()

	def post_and_parse_resp(self, payload):
		'''
		Do POST and parse response to pd.DataFrame
		'''
		resp = self.app.post(ENTRY_POINT + '/search', 
			data=json.dumps(payload),
			content_type = 'application/json')
		# Parse the repsonse data to pandas.DataFrame
		signatures = json.loads(resp.data.decode())
		signatures = pd.DataFrame.from_records(signatures).set_index('id')
		return resp, signatures

	def test_query(self):
		resp, signatures = self.post_and_parse_resp(self.payload)
		# Check response code and mimetype
		self.assertEquals(resp.status_code, 200)
		self.assertEquals(resp.mimetype, 'application/json')
		
		# Check all signed jaccard < 0 for opposite
		self.assertTrue(all(signatures['signed_jaccard'] < 0))

	def test_db_version(self):
		payload1 = self.payload.copy()
		payload2 = self.payload.copy()

		payload1['db_version'] = ['v1.0']
		payload2['db_version'] = ['v1.0','DM']

		resp, signatures = self.post_and_parse_resp(self.payload)
		resp1, signatures1 = self.post_and_parse_resp(payload1)
		resp2, signatures2 = self.post_and_parse_resp(payload2)
		self.assertEquals(resp2.status_code, 200)
		self.assertEquals(resp2.mimetype, 'application/json')
		# Make sure db_version = 'v1.0' is equivalent to ['v1.0']
		self.assertEquals(signatures1.shape[0], signatures.shape[0])

		# Make sure db_version = ['v1.0','DM'] has more results than v1.0
		self.assertTrue(signatures2.shape[0] > signatures.shape[0])


class TestRetrieveUsingId(unittest.TestCase):
	signature_id = 'gene:27'
	def setUp(self):
		self.app = app.test_client()

	def test_retrieve(self):
		resp = self.app.get(ENTRY_POINT + '/api?id=%s' % self.signature_id)

		self.assertEquals(resp.status_code, 200)
		self.assertEquals(resp.mimetype, 'application/json')

		# Parse the response object
		signature = json.loads(resp.data.decode())

		self.assertEquals(signature['id'], self.signature_id)
		self.assertEquals(type(signature['ctrl_ids']), list)
		self.assertEquals(type(signature['pert_ids']), list)



if __name__ == '__main__':
    unittest.main()
