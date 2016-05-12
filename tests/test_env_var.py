import os
import unittest

class TestEnvVar(unittest.TestCase):
	'''
	Test if environmental variable `CONFIG_OBJ` is set
	'''
	def test_env_var(self):
		self.assertEquals(os.environ['CONFIG_OBJ'], 'config.TestingConfig')

