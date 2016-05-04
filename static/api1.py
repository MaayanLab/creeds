import json, requests
from pprint import pprint

CREEDS_URL = 'http://amp.pharm.mssm.edu/CREEDS/'
response = requests.get(CREEDS_URL + 'search', params={'q':'STAT3'})
if response.status_code == 200:
	pprint(response.json())
	json.dump(response.json(), open('api1_result.json', 'wb'), indent=4)
