## the class for GEO entry
import urllib2
import json
import numpy as np
from scipy.stats import fisher_exact
from Bio import Entrez
Entrez.email = 'wangzc921@gmail.com'

class GEOentry(object):
	"""docstring for GEOentry"""
	def __init__(self, geo_id, ctrls, perts, gene, pert_type, platform, organism, cell, curator):
		ctrls.sort()
		perts.sort()
		self.geo_id = geo_id
		self.ctrls = ctrls
		self.perts = perts
		self.gene = gene # disease name for dz
		self.pert_type = pert_type # disease id for dz
		self.platform = platform
		self.organism = organism
		self.cell = cell
		self.curator = curator
		## more fields
		self.time = None
		self.uid = None
		self.chdir = None 
		self.conversion_pct = None


	def get_url(self, base_url="http://amp.pharm.mssm.edu/g2e/full?"):
		url = base_url
		## required args
		url += 'accession=%s&'%self.geo_id
		url += 'platform=%s&'%self.platform
		url += 'control=%s&'%'-'.join(self.ctrls)
		url += 'experimental=%s'%'-'.join(self.perts)
		## optional args
		# url += 'organism=%s&'%self.organism.lower().replace(' ', '+')
		# url += 'cell=%s&'%self.cell
		# url += 'perturbation=%s&'%self.pert_type.replace(' ', '+')
		# url += 'gene=%s'%self.gene
		url += '&cutoff=None' ## new args to retrieve the full chdir
		return url

	def get_lists(self):
		url = self.get_url()
		response = urllib2.urlopen(url)
		data = json.load(response)
		print 'API status:', data['status']
		dn_genes = data['down_genes'].keys()
		dn_chdir = map(float, data['down_genes'].values())
		srt_idx = np.argsort(dn_chdir)[::-1]
		dn_genes = np.array(dn_genes)[srt_idx]

		up_genes = data['up_genes'].keys()
		up_chdir = map(float, data['up_genes'].values())
		srt_idx = np.argsort(up_chdir)[::-1]
		up_genes = np.array(up_genes)[srt_idx]

		self.dn_genes = dn_genes
		self.up_genes = up_genes
		time.sleep(1)

	def get_json(self):
		url = self.get_url()
		response = urllib2.urlopen(url)
		data = json.load(response)
		print 'API status:', data['status']
		data['geo_id'] = self.geo_id
		data['ctrls'] = self.ctrls
		data['perts'] = self.perts
		data['gene'] = self.gene
		data['pert_type'] = self.pert_type
		data['platform'] = self.platform
		data['organism'] = self.organism
		data['cell'] = self.cell
		data['curator'] = self.curator
		return data


	def get_sql(self, table='geo2enrichr'):
		values = ','.join(["'%s'"%s for s in [self.geo_id, 
			','.join(self.ctrls), ','.join(self.perts), 
			self.gene, self.pert_type, self.organism, self.cell,
			','.join(self.up_genes), ','.join(self.dn_genes),
			self.curator] ] )
		sql = "INSERT INTO %s VALUES (%s, NOW(), NULL)" %(table, values)
		return sql

	def __eq__(self, entry):
		# to assert whether two entries are identical
		if self.ctrls == entry.ctrls and self.perts == entry.perts:
			return True
		else:
			return False

	def get_lists_cutoff(self, cutoff): ## get up/dn gene lists from chdir by applying rank cutoff
		if self.chdir is None:
			raise ValueError('chdir is not set!')
		else:
			chdir_values = np.array(self.chdir.values())
			genes = np.array(self.chdir.keys())
			srt_idx = np.argsort(chdir_values)
			genes = genes[srt_idx]
			self.dn_genes = genes[0: cutoff]
			self.up_genes = genes[-cutoff:]
			self.chdir = None # free the memory


def json2entry(fn):
	# retrieve entry from a json file
	json_data = json.load(open(fn,'r'))
	uid = int(fn.split('.')[0])
	entry = GEOentry(json_data['geo_id'], json_data['ctrls'], json_data['perts'], json_data['gene'], json_data['pert_type'], json_data['platform'], json_data['organism'], json_data['cell'], json_data['curator'])
	entry.uid = uid
	entry.conversion_pct = float(json_data['conversion_pct'])
	entry.time = json_data['time']
	# entry.chdir = dict( json_data['up_genes'].items() + json_data['dn_genes'].items() )
	chdir = {}
	for key, val in json_data['up_genes'].items() + json_data['down_genes'].items():
		chdir[key] = float(val)
	entry.chdir = chdir
	return entry

## algorithms
def fisher_exact_test(l1, l2, universe=22000):
	s1, s2 = set(l1), set(l2)
	overlap = len(s1 & s2)
	_, p = fisher_exact([[overlap, len(s1)], [len(s2), universe]])
	return p

def jaccard(l1,l2):
	s1, s2 = set(l1), set(l2)
	up = len(s1 & s2)
	dn = len(s1 | s2)
	return float(up)/dn

def signed_jaccard(e1, e2):
	## signed jaccard index 
	j1 = jaccard(e1.up_genes, e2.up_genes)
	j2 = jaccard(e1.dn_genes, e2.dn_genes)
	j3 = jaccard(e1.dn_genes, e2.up_genes)
	j4 = jaccard(e1.up_genes, e2.dn_genes)
	return (j1 + j2 - j3 - j4) / 2

# def cosine_dist(e1, e2):
# 	return

## utils
def geo_id2platform(geo_id):
	if geo_id.startswith('GDS'):
		geo_id = geo_id[3:]
		handle = Entrez.esummary(db='gds', id=geo_id)
		record = Entrez.read(handle)
		platform = 'GPL'+record[0]['GPL']
	else:
		handle = Entrez.esearch(db='gds', term='%s[GEO Accession]'%geo_id)
		records = Entrez.read(handle)
		for uid in records['IdList']:
			rec= Entrez.read(Entrez.esummary(db='gds', id=uid))[0]
			if 'GSE'+rec['GSE'] == geo_id:
				platform = 'GPL'+rec['GPL']
				break
	return platform
