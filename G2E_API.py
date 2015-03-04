## process the metadata of Kevin and Andrew
## and use GEO2Enrichr API to get gene lists
## created on 2/11/2015
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

conn = MySQLdb.connect(host='localhost',user='root', passwd='',db='maaya0_crowdsourcing')
cur = conn.cursor()

## 1. parse Andrew's meta file
# error_log = open('andrew5.log','w')
# # with open ('kinase_accession_list_20141112.txt') as f:
# with open ('andrew4.log') as f:
# 	for line in f:
# 		sl = line.strip().split('\t')
# 		geo_id = sl[0]
# 		gene, pert_type, uid = sl[1].split('_')
# 		ctrls = sl[2].split(',')
# 		perts = sl[3].split(',')
# 		platform = sl[4]
# 		organism = sl[5]
# 		cell = sl[6]

# 		try:
# 			entry = GEOentry(geo_id, ctrls, perts, gene, pert_type, platform, organism, cell, 'Andrew')
# 			print uid, entry.get_url()
# 			sql = entry.get_sql()
# 			try:
# 				cur.execute(sql)
# 				print 'success'
# 				conn.commit()
# 			except Exception, e:
# 				conn.rollback()
# 		except:
# 			print uid, 'error'
# 			error_log.write(line)


## 2. parse Kevin's meta file
# error_log = open('kevin5.log', 'w')
# # with open ('TFLOF_meta.txt') as f:
# with open ('kevin4.log') as f: 
# 	i = 0
# 	for line in f:
# 		i += 1
# 		uid = i
# 		sl = line.strip().split('\t')
# 		geo_id = sl[0]
# 		gene = sl[1]
# 		pert_type = sl[2]
		
# 		ctrls = sl[6].split(',')
# 		perts = sl[7].split(',')
# 		platform = sl[3]
# 		organism = sl[4]
# 		cell = sl[5]

# 		try:
# 			entry = GEOentry(geo_id, ctrls, perts, gene, pert_type, platform, organism, cell, 'Kevin')
# 			print uid, entry.get_url()
# 			sql = entry.get_sql()
# 			try:
# 				cur.execute(sql)
# 				print 'success'
# 				conn.commit()
# 			except Exception, e:
# 				conn.rollback()
# 		except:
# 			print uid, 'error'
# 			error_log.write(line)

## 3. parse entries from other curators in the database 
d_gds_gse = dict(zip(file2list('GDS_GSE_conversion.txt',0), file2list('GDS_GSE_conversion.txt', 1) ))
# query0 = """SELECT geo_id FROM geo2enrichr"""  # get all platform info
# cur.execute(query0)

# d_gse_platform = {}
# with open ('GSE_platform.txt', 'w') as out:
# 	for row in cur:
# 		geo_id = row[0]		
# 		if geo_id in d_gds_gse:
# 			geo_id = d_gds_gse[geo_id]
# 		if geo_id not in d_gse_platform:
# 			try:
# 				platform = geo_id2platform(geo_id)
# 			except Exception:
# 				platform = 'error'
# 				pass
# 			d_gse_platform[geo_id] = platform
# 	for gse in d_gse_platform:
# 		out.write(gse + '\t' + d_gse_platform[gse] + '\n')

d_gse_platform = dict(zip(file2list('GSE_platform.txt',0), file2list('GSE_platform.txt', 1) ))
'''
query1 = """SELECT * FROM geo2enrichr WHERE curator='Andrew' OR curator='Kevin'"""
cur.execute(query1)
valid_entries = {}

i = 0
for row in cur:
	i += 1
	geo_id, ctrls, perts, gene, pert_type, organism, cell = row[0:7]
	curator = row[-3]
	up_genes, dn_genes = row[7:9]
	uid = row[-1]

	if geo_id in d_gds_gse:
		geo_id = d_gds_gse[geo_id]
	platform = d_gse_platform[geo_id]
	entry = GEOentry(geo_id, ctrls.split(','), perts.split(','), gene, pert_type, platform, organism, cell, curator)
	valid_entries[uid] = entry
	if i % 200 == 0: print i

print "Andrew and Kevin's entries:", len(valid_entries)

query2 = """SELECT * FROM geo2enrichr WHERE curator!='Andrew' AND curator!='Kevin'"""
cur.execute(query2)

for row in cur:
	i += 1
	geo_id, ctrls, perts, gene, pert_type, organism, cell = [item.strip() for item in row[0:7]]
	geo_id = geo_id.strip().split('[')[0]
	uid = row[-1]
	curator = row[-3]

	if geo_id in d_gds_gse:
		geo_id = d_gds_gse[geo_id]
	if geo_id in d_gse_platform:
		platform = d_gse_platform[geo_id]
	else:
		platform = 'error'
	entry = GEOentry(geo_id, ctrls.split(','), perts.split(','), gene, pert_type, platform, organism, cell, curator)
	bools = [valid_e == entry for valid_e in valid_entries.values()]
	if sum(bools) == 0: # all False
		valid_entries[uid] = entry
	if i % 200 == 0: print i

print "All unique entries:", len(valid_entries)
'''
## use G2E API to get list and dump them to json
# os.chdir('D:\Zichen_Projects\microtask_GEO\jsons')
# exist_entries = os.listdir('D:\Zichen_Projects\microtask_GEO\jsons')
# for uid, entry in valid_entries.items():
# 	fn = str(uid) + '.json'
# 	if fn not in exist_entries:
# 		url = entry.get_url()
# 		print url
# 		try:
# 			json_data = entry.get_json()
# 			json.dump(json_data, open(fn, 'w'))
# 		except Exception:
# 			print 'error'
# 			pass


# pickle.dump(valid_entries, open('valid_entries.p','w'))


# valid_entries = pickle.load(open('valid_entries.p','r')) ## uniq entries with platform info
# valid_entries = pickle.load(open('entries_error_in_API.p','r')) ## uniq entries with platform info
# print 'number of valid_entries:', len(valid_entries)
# entries_error_in_API = []

# ## re-run G2E API to get lists for entries from the crowd and put in to a new table `geo2enrichr_cleaned`
# for e in valid_entries:
# 	if e.platform == 'error':
# 		print e.geo_id, e.ctrls, e.perts
# 	else:
# 		if e.curator in ['Andrew', 'Kevin']:
# 			sql = e.get_sql(table='geo2enrichr_cleaned')
# 			cur.execute(sql)
# 			conn.commit()
# 		else: ## entries from the crowd
# 			try:
# 				print e.get_url()
# 				e.get_lists()
# 				sql = e.get_sql(table='geo2enrichr_cleaned')
# 				try:
# 					cur.execute(sql)
# 					print 'success'
# 					conn.commit()
# 				except Exception, e:
# 					conn.rollback()
# 			except:
# 				entries_error_in_API.append(e)

# # error_log.close()
# pickle.dump(entries_error_in_API, open('entries_error_in_API2.p','w'))


## use G2E API to process nosology
# query3 = """SELECT * FROM geo2enrichr_dz"""
# cur.execute(query3)
# valid_entries = {}
# for row in cur:
# 	geo_id, ctrls, perts, gene, pert_type, organism, cell = [item.strip() for item in row[0:7]]
# 	uid = row[-1]
# 	curator = row[-3]

# 	if geo_id in d_gds_gse:
# 		geo_id = d_gds_gse[geo_id]
# 	if geo_id in d_gse_platform:
# 		platform = d_gse_platform[geo_id]
# 		entry = GEOentry(geo_id, ctrls.split(','), perts.split(','), gene, pert_type, platform, organism, cell, curator)
# 		bools = [valid_e == entry for valid_e in valid_entries.values()]
# 		if sum(bools) == 0: # all False
# 			valid_entries[uid] = entry


# pickle.dump(valid_entries, open('valid_entries_dz.p', 'w'))
valid_entries = pickle.load(open('valid_entries_dz.p','r'))
print "All unique entries:", len(valid_entries)

os.chdir('D:\Zichen_Projects\microtask_GEO\jsons_dz')
exist_entries = os.listdir('D:\Zichen_Projects\microtask_GEO\jsons_dz')
for uid, entry in valid_entries.items():
	fn = str(uid) + '.json'
	if fn not in exist_entries:
		url = entry.get_url()
		print url
		try:
			json_data = entry.get_json()
			json.dump(json_data, open(fn, 'w'))
		except Exception:
			print 'error', entry.platform
			pass

conn.close()
