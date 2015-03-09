## run the API of G2E of branch `zichen` to re-analyze the GEO entries
## created on 3/9/2015
import os, sys
import numpy as np
import time
import urllib2
import json
import MySQLdb
import cPickle as pickle
from pprint import pprint

from GEOentry import GEOentry, geo_id2platform

from Bio import Entrez
Entrez.email = 'wangzc921@gmail.com'
sys.path.append('C:\Users\Zichen\Documents\\bitbucket\maayanlab_utils')
from fileIO import file2list, mysqlTable2dict

##
DATADIR = 'output/annot_jsons/'

## get dicts from mysql
d_uid_hsgene = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_genes', 0, 1)
d_uid_mmgene = mysqlTable2dict('maaya0_crowdsourcing', 'cleaned_genes', 0, 2)
d_gds_gse = mysqlTable2dict('maaya0_crowdsourcing', 'gds_gse', 0, 1)
d_gse_platform = mysqlTable2dict('maaya0_crowdsourcing', 'geo_platforms', 0, 1)


conn = MySQLdb.connect(host='localhost',user='root', passwd='',db='maaya0_crowdsourcing')
cur = conn.cursor()
cur_insert = conn.cursor()

query = """SELECT * FROM geo2enrichr"""
cur.execute(query)

i = 0
valid_entries = {}

for row in cur:
	i += 1
	# print row
	geo_id, ctrls, perts, gene, pert_type, organism, cell = [item.strip() for item in row[0:7]]
	geo_id = geo_id.strip().split('[')[0]
	uid = row[-1]
	curator = row[-3]
	## get cleaned gene symbol
	if organism in ['mouse', 'Mus musculus']:
		if d_uid_mmgene[uid] != 'NULL':
			gene = d_uid_mmgene[uid]
		else:
			gene = d_uid_hsgene[uid]
	elif organism in ['human', 'Homo sapiens']:
		if d_uid_hsgene[uid] != 'NULL':
			gene = d_uid_hsgene[uid]
		else:
			gene = d_uid_mmgene[uid]
	## get platform
	if geo_id in d_gds_gse:
		geo_id = d_gds_gse[geo_id]
	if geo_id in d_gse_platform:
		platform = d_gse_platform[geo_id]
	else:
		try:
			platform = geo_id2platform(geo_id)
			cur_insert.execute("INSERT INTO geo_platforms VALUES (%s, %s)", (geo_id, platform))
		except:
			platform = 'error'
			print uid, row
			pass
	if platform != 'error' and gene != 'NULL':
		entry = GEOentry(geo_id, filter(None,ctrls.split(',')), filter(None,perts.split(',')), gene, pert_type, platform, organism, cell, curator)
		bools = [valid_e == entry for valid_e in valid_entries.values()]
		if sum(bools) == 0: # all False
			valid_entries[uid] = entry
	if i % 200 == 0: print i

conn.commit()
print "All unique entries:", len(valid_entries)
conn.close()

## use G2E API to get list and dump them to json
os.chdir(DATADIR)
error_log = open('errors.log', 'w')
exist_entries = os.listdir(os.getcwd())
for uid, entry in valid_entries.items():
	fn = str(uid) + '.json'
	if fn not in exist_entries:
		url = entry.get_url('http://127.0.0.1:8083/g2e/full?')
		print uid, url
		try:
			json_data = entry.get_json()
			json.dump(json_data, open(fn, 'w'))
		except Exception, e:
			print uid, 'error', e
			error_log.write(str(uid) + '\t' + e.message + '\t' + url + '\n')
			pass
