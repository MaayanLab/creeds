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
	# list of parameters for DBSignatureCollection globals
	DBSC_PARAMS = [
		(
			{'$and':[ 
				{'chdir_sva_exp2': {'$exists': True}}, 
				{'version': '1.0'},
				{"incorrect": {"$ne": True}},
			]},
			'v1.0'),

		(
			{'$and':[
				{'chdir_sva_exp2': {'$exists': True}}, 
				{'version': '1.2'},
			]},
			'DM'
			)
	]


class ProductionConfig(Config):
	DATABASE_URI = 'mongodb://146.203.54.131:27017/'
	HOST = '0.0.0.0'


class DevelopmentConfig(Config):
	DEBUG = True


class TestingConfig(DevelopmentConfig):
	DBSC_PARAMS = [
		(
			{'$and':[ 
				{'chdir_sva_exp2': {'$exists': True}}, 
				{'version': '1.0'},
				{"incorrect": {"$ne": True}},
				{'id': {'$in': ['gene:27', 'gene:3046', 'gene:2981', 'gene:1829']}}
			]},
			'v1.0'),

		(
			{'$and':[
				{'chdir_sva_exp2': {'$exists': True}}, 
				{'version': '1.2'},
				{'id': {'$in': ['drug:DM0', 'drug:DM1', 'drug:DM10']}}
			]},
			'DM'
			)
	]

		