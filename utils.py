'''
Utils for computing similarity and etc.
'''
def get_matrix(uids, genes, na_val=0):
	## retrieve a matrix based on uids of signatures and genes
	mat = np.zeros((len(genes), len(uids)))

	for j, uid in enumerate(uids):
		sig = DBSignature(uid, projection=PROJECTION_EXCLUDE)
		vals = sig.get_gene_vals(genes, na_val=na_val)
		mat[:, j] = vals
	return mat


def make_autocomplete(db_collections):
	'''
	Make the object required for autocomplete 
	`db_collections`: a list of DBSignatureCollection instances
	'''
	d_cat_names = {}
	for uid, sig in [dbc.items() for dbc in db_collections]:
		if sig.has_chdir():
			cat = uid.split(':')[0]
			if cat not in d_cat_names:
				d_cat_names[cat] = []
			name = sig.name
			d_cat_names[cat].append(name)

	for cat, names in d_cat_names.items():
		d_cat_names[cat] = list(set(names))

	return d_cat_names


