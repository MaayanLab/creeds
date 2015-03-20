## script to parse file from diseaes-ontology 
## 'HumanDO.obo'
## created on 3/20/2015

import os, sys
sys.path.append('C:\\Users\\Zichen\\Documents\\bitbucket\\maayanlab_utils')
import obo_parser


DO_FILE = 'D:\Zichen_Projects\microtask_GEO\DO\HumanDO.obo'
ontology_id_terms, name2id, G = obo_parser.obo2network(DO_FILE)
print len(ontology_id_terms), len(name2id), G.number_of_edges(), G.number_of_nodes()
# print ontology_id_terms['DOID:0001816'].xrefs
# print ontology_id_terms['DOID:0001816'].synonym
# print ontology_id_terms['DOID:0001816'].name
# print ontology_id_terms['DOID:0001816'].parents
# print ontology_id_terms['DOID:0001816'].children

d_doid_umls_ids = {}
d_umls_doid = {}

for doid, term in ontology_id_terms.items():
	umls_ids = [xref for xref in term.xrefs if xref.startswith('UMLS_CUI:')]
	d_doid_umls_ids[doid] = umls_ids
	for umls_id in umls_ids:
		d_umls_doid[umls_id] = doid

# output the dictionary
with open ('umls_id2_DOID.info', 'w') as out:
	for umls_id, doid in d_umls_doid.items():
		out.write(umls_id+'\t'+doid+'\n')

