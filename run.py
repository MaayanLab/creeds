import os, sys


if __name__ == '__main__':
	if len(sys.argv) > 1:
		mode = sys.argv[1]
	else:
		mode = 'dev'

	if mode == 'dev':
		os.environ['CONFIG_OBJ'] = 'config.DevelopmentConfig'
	elif mode == 'test':
		os.environ['CONFIG_OBJ'] = 'config.TestingConfig'
	elif mode == 'speed_test':
		os.environ['CONFIG_OBJ'] = 'config.SpeedTestingConfig'
	else:
		os.environ['CONFIG_OBJ'] = 'config.ProductionConfig'

	from creeds import app

	host = app.config['HOST']
	port = app.config['PORT']

	app.run(host=host, port=port)

