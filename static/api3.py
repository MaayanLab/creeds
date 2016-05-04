import json, requests
from pprint import pprint

CREEDS_URL = 'http://amp.pharm.mssm.edu/CREEDS/'
response = requests.get(CREEDS_URL + 'api', params={'id':'gene:84'})
if response.status_code == 200:
	pprint(response.json())
	json.dump(response.json(), open('api3_result.json', 'wb'), indent=4)
