## Parse Axel's GEO metadata

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

os.chdir('D:\Zichen_Projects\drug_se_prediction\DrugBank')
d_name_dbid = dict(zip( file2list('all_DrugBank_pert_id.info', 1), file2list('all_DrugBank_pert_id.info', 0) ))

os.chdir('D:\Zichen_Projects\microtask_GEO')

conn = MySQLdb.connect(host='localhost',user='root', passwd='',db='maaya0_crowdsourcing')
cur = conn.cursor()

with open ('GDSfinalOutput.txt') as f:
	for line in f:
		sl = line.strip().split('\t')
		geo_id = sl[0].upper()
		drug = sl[1].lower()
		if drug.upper() in d_name_dbid:
			drug_id = d_name_dbid[drug.upper()]
		else:
			drug_id = 'NA'
		ctrls = [item.strip().upper() for item in sl[2].split(',')]
		perts = [item.strip().upper() for item in sl[3].split(',')]
		ctrls.sort()
		perts.sort()
		if sl[4] == '[]':
			cell = 'NA'
		else:
			cell = sl[4]
		organism = sl[5]
		platform = sl[6]
		values = ','.join( ["'%s'"%item for item in [geo_id, ','.join(ctrls), ','.join(perts), drug, drug_id, organism, cell,'NA','NA', 'Axel'] ])
		sql = """INSERT INTO geo2enrichr_drug VALUES (%s, NOW(), NULL)"""%values
		print sql
		try:
			cur.execute(sql)
			# print sql
			conn.commit()
		except Exception as e:
			print e
			conn.rollback()
conn.close()
