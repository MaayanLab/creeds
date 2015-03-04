## get GEO platform info and put into MySQL database

import os, sys
import numpy as np
import time
import urllib2
import json
import MySQLdb
import cPickle as pickle
from pprint import pprint

from GEOentry import GEOentry

from Bio import Entrez
Entrez.email = 'wangzc921@gmail.com'
sys.path.append('C:\Users\Zichen\Documents\\bitbucket\maayanlab_utils')
from fileIO import file2list

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


os.chdir('D:\Zichen_Projects\microtask_GEO')
d_gds_gse = dict(zip(file2list('GDS_GSE_conversion.txt',0), file2list('GDS_GSE_conversion.txt', 1) ))
d_gse_platform = dict(zip(file2list('GSE_platform.txt',0), file2list('GSE_platform.txt', 1) ))

conn = MySQLdb.connect(host='localhost',user='root', passwd='',db='maaya0_crowdsourcing')
cur = conn.cursor()

query3 = """SELECT * FROM geo2enrichr_dz"""
cur.execute(query3)

out = open('GSE_platform.txt', 'a+')

for row in cur:
	geo_id, ctrls, perts, gene, pert_type, organism, cell = [item.strip() for item in row[0:7]]
	uid = row[-1]
	curator = row[-3]

	if geo_id in d_gds_gse:
		geo_id = d_gds_gse[geo_id]
	if geo_id in d_gse_platform:
		platform = d_gse_platform[geo_id]
	else:
		try:
			platform = geo_id2platform(geo_id)
			print geo_id, platform
			out.write(geo_id +"\t"+ platform+' \n')
		except:
			out.write(geo_id +"\t"+ 'error'+' \n')

out.close()