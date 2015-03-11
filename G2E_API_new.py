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
		ctrls = filter(None,ctrls.split(','))
		perts = filter(None,perts.split(','))
		if len(ctrls) == 1 or len(perts) == 1: 
			if ';' in ctrls[0] or ';' in perts[0]:
				delim = ';'
			else:
				delim = ' '
			ctrls = [s.strip() for s in ctrls[0].split(delim)]
			perts = [s.strip() for s in perts[0].split(delim)]

		if len(ctrls) != 1 and len(perts) != 1: # should have replicates
			if len(set(ctrls) & set(perts)) == 0: # GSMs are not allowed to be in both ctrls and perts
				entry = GEOentry(geo_id, ctrls, perts, gene, pert_type, platform, organism, cell, curator)
				bools = [valid_e == entry for valid_e in valid_entries.values()]
				if sum(bools) == 0: # all False
					valid_entries[uid] = entry
			else:
				print uid
	# if i % 200 == 0: print i

print "All unique entries:", len(valid_entries)
# os.chdir(DATADIR)
## run once to put all valid entries into table `valid_microtask1`
# no_gsm_in_soft_errors = map(int, file2list('errors.log', 0))
# for uid in valid_entries:
# 	if uid not in no_gsm_in_soft_errors:
# 		entry = valid_entries[uid]
# 		cur_insert.execute('INSERT INTO valid_microtask1 VALUES (%s, %s, %s, %s, %s)', 
# 			(uid, entry.geo_id, ','.join(entry.ctrls), ','.join(entry.perts), entry.platform )
# 			)

conn.commit()
conn.close()

## use G2E API to get list and dump them to json
BASE_URL = 'http://127.0.0.1:8083/g2e/full?'
os.chdir(DATADIR)
# error_log = open('errors.log', 'w')
# exist_entries = os.listdir(os.getcwd())
# for uid, entry in valid_entries.items():
# 	fn = str(uid) + '.json'
# 	if fn not in exist_entries:
# 		url = entry.get_url(BASE_URL)
# 		print uid, url
# 		try:
# 			json_data = entry.get_json(BASE_URL)
# 			json.dump(json_data, open(fn, 'w'))
# 		except Exception, e:
# 			print uid, 'error', e
# 			error_log.write(str(uid) + '\t' + e.message + '\t' + url + '\n')
# 			pass
# error_log.close()

## re-run unsupported platforms that don't have annot files 3/11/2015
error_log = open('errors3.log', 'w')
with open ('no_genes_errors.log') as f:
	for line in f:
		sl = line.strip().split('\t')
		uid = int(sl[0])
		if sl[1] == 'error' and uid in valid_entries:
			fn = str(uid) + '.json'
			entry = valid_entries[uid]
			url = entry.get_url(BASE_URL)
			try:
				json_data = entry.get_json(BASE_URL)
				json.dump(json_data, open(fn, 'w'))
			except Exception, e:
				print uid, 'error', e
				error_log.write(str(uid) + '\t' + e.message + '\t' + url + '\n')
				pass
error_log.close()
