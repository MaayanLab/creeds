## to do HC and output a json compatible with Nick's visualization

import numpy as np
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
from scipy.stats import zscore


def clustergram(data, rids, cids,
	row_linkage='average', col_linkage='average', 
	row_pdist='euclidean', col_pdist='euclidean',
	standardize=3, log=False):
	## preprocess data
	if log:
		data = np.log2(data + 1.0)

	if standardize == 1: # Standardize along the columns of data
		data = zscore(data, axis=0)
	elif standardize == 2: # Standardize along the rows of data
		data = zscore(data, axis=1)

	## perform hierarchical clustering for rows and cols
	## compute pdist for rows:
	d1 = dist.pdist(data, metric=row_pdist)
	D1 = dist.squareform(d1)
	Y1 = sch.linkage(D1, method=row_linkage, metric=row_pdist)
	Z1 = sch.dendrogram(Y1, orientation='right')
	idx1 = Z1['leaves']

	## compute pdist for cols
	d2 = dist.pdist(data.T, metric=col_pdist)
	D2 = dist.squareform(d2)
	Y2 = sch.linkage(D2, method=col_linkage, metric=col_pdist)
	Z2 = sch.dendrogram(Y2)
	idx2 = Z2['leaves']

	row_nodes = []
	for idx, rid in zip(idx1, rids):
		row_nodes.append({'sort': idx, 'name': rid})

	col_nodes = []
	for idx, cid in zip(idx2, cids):
		col_nodes.append({'sort': idx, 'name': cid})	

	links = []
	for i in range(len(rids)):
		for j in range(len(cids)):
			links.append({'source': i, 'target': j, 'value': data[i,j]})
	
	json_data = {
		'row_nodes':row_nodes,
		'col_nodes':col_nodes,
		'links': links
				}
	return json_data


## to test
# np.random.seed(1)
# data=np.random.randn(3,6)
# # data[0,0] = np.nan
# # print data
# cids=['a','b','c','d','e','f']
# rids=['1','2','3']

# json_data = clustergram(data, rids, cids)
# # print json_data
# json.dump(json_data, open('data/test_clustergram.json', 'wb'))


