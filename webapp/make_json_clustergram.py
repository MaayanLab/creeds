import os, sys, json
from clustergram import clustergram
from orm_utils import *
sys.path.append('../../maayanlab_utils')
from fileIO import read_df




def find_name_for_id(uid):
	projection = {'_id':False, 'limma':False, 'fold_changes':False, 'chdir': False}
	dbs = DBSignature(uid,projection)
	return dbs.name + '|' + dbs.meta['geo_id']

# mat, ids, _ = read_df('../signed_jaccard_subset_unique_entries_831x831.txt')
# mat, ids, _ = read_df('../signed_jaccard_subset_unique_entries_519x519.txt')
mat, ids, _ = read_df('../signed_jaccard_subset_unique_entries_259x259.txt')

names = map(find_name_for_id, ids)

json_data = clustergram(mat, names, names,
	row_linkage='average', col_linkage='average',
	row_pdist='cosine', col_pdist='cosine')

json.dump(json_data, open('/Library/WebServer/Documents/d3_clustergram/signed_jaccard_subset_clustergram_259.json', 'wb'))	


