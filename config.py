'''
Config objects.
'''

class Config(object):
	"""Default configs"""
	DEBUG = False
	SEND_FILE_MAX_AGE_DEFAULT = 6
	DATABASE_URI = 'mongodb://127.0.0.1:27017/'
	HOST = '127.0.0.1'
	PORT = 5000
	# list of kwargs for DBSignatureCollection globals
	DBSC_PARAMS = [
		{
			'filter_': {'$and':[ 
				{'chdir_sva_exp2': {'$exists': True}}, 
				{'version': '1.0'},
				{"incorrect": {"$ne": True}},
			]},
			'name': 'v1.0', 'name_prefix': 'Manual'},

		{
			'filter_': {'$and':[
				{'chdir_sva_exp2': {'$exists': True}}, 
				{'version': '1.2'},
			]},
			'name': 'DM', 'name_prefix': 'DrugMatrix'
			},
		{
			'filter_': {'$and':[
				{'chdir_sva_exp2': {'$exists': True}},
				{'version': '2.0'},
			]},
			'name': 'p1.0', 'name_prefix': 'Automatic'
			}
	]
	# whether to generate files for downloading from DBSC instances
	MAKE_DOWNLOAD_FILES = True


class ProductionConfig(Config):
	DATABASE_URI = 'mongodb://146.203.54.131:27017/'
	HOST = '0.0.0.0'


class DevelopmentConfig(Config):
	DEBUG = True


class TestingConfig(DevelopmentConfig):
	DBSC_PARAMS = [
		{
			'filter_': {'$and':[ 
				{'chdir_sva_exp2': {'$exists': True}}, 
				{'version': '1.0'},
				{"incorrect": {"$ne": True}},
				{'id': {'$in': ['gene:27', 'gene:3046', 'gene:2981', 'gene:1829']}}
			]},
			'name': 'v1.0', 'name_prefix': 'Manual'},

		{
			'filter_': {'$and':[
				{'chdir_sva_exp2': {'$exists': True}}, 
				{'version': '1.2'},
				{'id': {'$in': ['drug:DM0', 'drug:DM1', 'drug:DM10', 'drug:DM11', 'drug:DM12']}}
			]},
			'name': 'DM', 'name_prefix': 'DrugMatrix'
			},
		{
			'filter_': {'$and':[
				{'chdir_sva_exp2': {'$exists': True}},
				{'version': '2.0'},
				{'id': {'$in': ['gene:P9030', 'dz:P1814', 'drug:P2341']}}
			]},
			'name': 'p1.0', 'name_prefix': 'Automatic'
			}
	]



class SpeedTestingConfig(DevelopmentConfig):
	'''
	Used for speed test.
	'''
	DBSC_PARAMS = [
		{
			'filter_': {'$and':[ 
				{'chdir_sva_exp2': {'$exists': True}}, 
				{'version': '1.0'},
				{"incorrect": {"$ne": True}},
				{'id': {'$in': ['gene:27', 'gene:3046', 'gene:2981', 'gene:1829']}}
			]},
			'name': 'v1.0', 'name_prefix': 'Manual',
			'limit': 300
			},

		{
			'filter_': {'$and':[
				{'chdir_sva_exp2': {'$exists': True}}, 
				{'version': '1.2'},
				{'id': {'$in': ['drug:DM0', 'drug:DM1', 'drug:DM10', 'drug:DM11', 'drug:DM12']}}
			]},
			'name': 'DM', 'name_prefix': 'DrugMatrix',
			'limit': 100
			},
		{
			'filter_': {'$and':[
				{'chdir_sva_exp2': {'$exists': True}},
				{'version': '2.0'},
				{'id': {'$in': ['gene:P9030', 'dz:P1814', 'drug:P2341']}}
			]},
			'name': 'p1.0', 'name_prefix': 'Automatic',
			'limit': 300
			}
	]

