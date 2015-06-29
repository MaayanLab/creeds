## the class for GEO entry
import urllib2
import json
import numpy as np
from scipy.stats import fisher_exact
from collections import OrderedDict
from Bio import Entrez
Entrez.email = 'wangzc921@gmail.com'

from gene_convertion import *

ORGANISMS = {
	'human': 'human', 'homo sapiens':'human',
	'mouse': 'mouse', 'mus musculus': 'mouse',
	'rat': 'rat', 'rattus norvegicus': 'rat',
} # control vocabularies for organisms

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
		self.organism = ORGANISMS[organism.lower()]
		self.cell = cell
		self.curator = curator
		## more fields
		self.time = None
		self.uid = None
		self.chdir = None 
		self.conversion_pct = None
		self.status = None
		self.message = None
		self.failed_to_download = None


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

	def get_json(self, base_url="http://amp.pharm.mssm.edu/g2e/full?"):
		url = self.get_url(base_url)
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

	def __len__(self):
		# the number of replicates
		return len(self.ctrls) + len(self.perts)

	def get_lists_cutoff(self, cutoff, to_human=False): ## get up/dn gene lists from chdir by applying rank cutoff
		if self.chdir is None:
			raise ValueError('chdir is not set!')
		else:
			chdir_values = np.array(self.chdir.values())
			genes = np.array(self.chdir.keys())
			self.chdir = None # free the memory
			srt_idx = np.argsort(chdir_values)
			genes = genes[srt_idx]
			if not to_human:
				self.dn_genes = clean_genes(genes[0: cutoff])
				self.up_genes = clean_genes(genes[-cutoff:])
			else:
				self.dn_genes = humanize(genes[0: cutoff])
				self.up_genes = humanize(genes[-cutoff:])

	def to_json_geneset(self, cutoff=None, fuzzy=False):
		# writen for dzs
		if self.chdir is None:
			raise ValueError('chdir is not set!')
		else:
			geneset = OrderedDict()
			if not fuzzy:
				self.get_lists_cutoff(cutoff, to_human=True)
				geneset['term'] = self.gene + '_' + self.geo_id
				geneset['desc'] = self.cell
				geneset['up'] = self.up_genes
				geneset['dn'] = self.dn_genes

			else: # fuzzy
				chdir_values = np.array(self.chdir.values())
				genes = np.array(self.chdir.keys())
				srt_idx = np.argsort(chdir_values)
				genes = humanize(genes[srt_idx])

				if cutoff is None:
					geneset['term'] = self.gene + '_' + self.geo_id
					geneset['desc'] = self.cell
					geneset['genes'] = genes
					geneset['vals'] = chdir_values[srt_idx].tolist()

				else:
					chdir_values = chdir_values[srt_idx]
					geneset['term'] = self.gene + '_' + self.geo_id
					geneset['desc'] = self.cell
					geneset['genes'] = genes[0:cutoff] + genes[-cutoff:]
					geneset['vals'] = chdir_values.tolist()[0:cutoff] + chdir_values.tolist()[-cutoff:]
			return geneset

	def to_full_chdir(self): ## for the PAEA shiny app
		geneset = OrderedDict()
		chdir_values = self.chdir.values()
		genes = self.chdir.keys()
		genes, chdir_values = clean_genes2(genes, chdir_values)
		geneset['term'] = self.gene + '_' + self.geo_id
		geneset['desc'] = self.cell
		geneset['genes'] = genes
		geneset['vals'] = chdir_values
		return geneset

	def to_json_geneset2(self, cutoff): 
		geneset = OrderedDict()
		chdir_values = self.chdir.values()
		genes = self.chdir.keys()
		genes, chdir_values = clean_genes2(genes, chdir_values)
		genes, chdir_values = np.array(genes), np.array(chdir_values)
		srt_idx = np.argsort(abs(chdir_values))[::-1]
		# geneset['term'] = self.uid
		geneset['term'] = '_'.join([self.gene, self.pert_type, self.organism, self.geo_id])
		geneset['genes'] = genes[srt_idx][0:cutoff].tolist()
		geneset['vals'] = chdir_values[srt_idx][0:cutoff].tolist()
		return geneset

	def to_sorted_gene_list(self, absolute=True):
		# return an array of cleaned genes sorted by chdir coefs
		# in descending order
		chdir_values = self.chdir.values()
		genes = self.chdir.keys()
		genes, chdir_values = clean_genes2(genes, chdir_values)
		genes, chdir_values = np.array(genes), np.array(chdir_values)
		if absolute:
			srt_idx = np.argsort(abs(chdir_values))[::-1]
		else:
			srt_idx = np.argsort(chdir_values)[::-1]
		return genes[srt_idx]

def json2entry(fn, meta_only=False):
	# retrieve entry from a json file
	json_data = json.load(open(fn,'r'))
	uid = int(fn.split('.')[0])
	entry = GEOentry(json_data['geo_id'], json_data['ctrls'], json_data['perts'], json_data['gene'], json_data['pert_type'], json_data['platform'], json_data['organism'], json_data['cell'], json_data['curator'])
	entry.uid = uid
	try:
		entry.conversion_pct = float(json_data['conversion_pct'])
	except:
		pass
	entry.time = json_data['time']
	if not meta_only:
		chdir = {}
		for key, val in json_data['up_genes'].items() + json_data['down_genes'].items():
			chdir[key] = float(val)
		entry.chdir = chdir
	else:
		if 'up_genes' in json_data:
			entry.chdir = len(json_data['up_genes']) + len(json_data['down_genes'])
		else:
			entry.chdir = 0
		entry.status = json_data['status']
		try:
			entry.message = json_data['message']
			entry.failed_to_download = json_data['failed_to_download']
		except:
			pass
	return entry

def json2genes(fn):
	# retrieve genes measured in an entry
	json_data = json.load(open(fn,'r'))
	uid = int(fn.split('.')[0])
	genes = json_data['up_genes'].keys() + json_data['down_genes'].keys()
	return genes

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
