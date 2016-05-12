import os, json
import unittest
from flask import (request, Response)
# Assumes env var for config is set
from creeds import app


class TestApiRequests(unittest.TestCase):
	ENTRY_POINT = '/creeds'
	def test_get_request_args(self):
		with app.test_request_context(self.ENTRY_POINT +'/api?id=gene:27'):
			self.assertEquals(request.path, self.ENTRY_POINT + '/api')
			self.assertEquals(request.args['id'], 'gene:27')



class TestApiResponseBase(unittest.TestCase):
	'''
	Base TestCase for testing response from API endpoints
	'''
	ENTRY_POINT = '/creeds'
	def setUp(self):
		self.app = app.test_client()


class TestStringSearch(TestApiResponseBase):
	'''
	Test the API of search the DB using string
	'''
	def test_search(self):
		# Send GET request and get Response object
		resp = self.app.get(self.ENTRY_POINT + '/search?q=TP53')

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


# class Test
