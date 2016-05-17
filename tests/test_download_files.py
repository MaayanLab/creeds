import os, json
import unittest
from StringIO import StringIO
import pandas as pd
from flask import (request, Response)
# Assumes env var for config is set
from creeds import app

ENTRY_POINT = '/creeds'

class TestDownloadFiles(unittest.TestCase):
	'''
	Test the files were able to be downloaded
	'''
	endpoint = ENTRY_POINT + '/download/'
	def setUp(self):
		self.dir = os.path.dirname(__file__)
		self.app = app.test_client()

	def download_file(self, filename):
		resp = self.app.get(self.endpoint + filename,
			follow_redirects=True)
		
		self.assertEquals(resp.status_code, 200)
		self.assertEquals(type(resp.data), str)
		return resp

	def test_downloaded_csv(self):
		resp = self.download_file('Single_gene_perturbations-v1.0.csv')
		self.assertEquals(resp.mimetype, 'text/csv')

		# Parse the response into pd.DataFrame
		resp_df = pd.read_csv(StringIO(resp.data))
		# Expect to have 4 signatures in this file
		self.assertEquals(resp_df.shape[0], 4)

	def test_downloaded_json(self):
		resp = self.download_file('Single_gene_perturbations-v1.0.json')
		self.assertEquals(resp.mimetype, 'application/json')
		# Parse
		resp = json.loads(resp.data)
		# Check data types of the parsed resp is list of dicts
		self.assertEquals(type(resp), list)
		self.assertEquals(len(resp), 4)
		self.assertEquals(type(resp[0]), dict)

	def test_downloaded_gmt(self):
		resp = self.download_file('Single_gene_perturbations-v1.0.gmt')
		self.assertEquals(resp.mimetype, 'application/octet-stream')
		# Parse
		d_gmt = {}
		for line in StringIO(resp.data):
			sl = line.strip().split('\t')
			key = '|'.join(sl[0:2])
			d_gmt[key] = sl[2:]

		self.assertEquals(len(d_gmt), 8)



if __name__ == '__main__':
	unittest.main()
