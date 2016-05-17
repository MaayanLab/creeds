## To test the speed for up/down gene set queries
import os
os.environ['CONFIG_OBJ'] = 'config.TestingConfig'
import numpy as np
import requests

from creeds.orm import GENE_SYMBOLS
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 
rcParams['font.sans-serif'] = 'Arial'
np.random.seed(42)

def generate_random_signature(size):
	## make a random signatures of with equal number of up and down genes
	all_genes = np.array(GENE_SYMBOLS.keys())
	rand_idx = np.random.randint(0, len(GENE_SYMBOLS), size=size * 2)
	rand_genes = all_genes[rand_idx]
	up_genes = rand_genes[:size].tolist()
	dn_genes = rand_genes[size:].tolist()
	return up_genes, dn_genes

def query_signatures(up_genes, dn_genes):
	## hit the POST API endpoint to query signatures with up/down genes
	# URL = 'http://amp.pharm.mssm.edu/CREEDS/search'
	URL = 'http://127.0.0.1:5000/CREEDS/search'
	res = requests.post(URL, 
		json={'up_genes': up_genes, 'dn_genes': dn_genes, 'direction':'similar'})
	assert res.status_code == 200
	return res

if __name__ == '__main__':
	import timeit

	sizes = [10, 50, 100, 200, 300, 500, 1000, 1500, 2000]
	times = np.zeros((len(sizes), 5), dtype=np.float32)
	for i, size in enumerate(sizes):

		setup = '''
from __main__ import query_signatures, generate_random_signature
up_genes, dn_genes = generate_random_signature(%s)
		''' % size
		t = timeit.Timer(stmt='query_signatures(up_genes, dn_genes)', 
			setup=setup).repeat(5, 5)	
		print size, t
		
		times[i] = t
	
	times /= 5
	plt.errorbar(sizes, times.mean(axis=1), yerr=times.std(axis=1), fmt='-o')
	plt.xlim([0, 2010])
	plt.xlabel('Size of gene set', fontsize=20)
	plt.ylabel('Time (seconds)', fontsize=20)
	plt.show()
