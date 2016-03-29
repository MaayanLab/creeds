import requests
BASE_URL = 'http://amp.pharm.mssm.edu/CREEDS/'

## 1. search by string
r = requests.get(BASE_URL + 'search', params={'q':'breast cancer'})
response = r.json()

print r.status_code
print response[0].keys()

## 2. search using gene sets
up_genes = ['KIAA0907','KDM5A','CDC25A','EGR1','GADD45B','RELB','TERF2IP','SMNDC1','TICAM1','NFKB2','RGS2','NCOA3','ICAM1','TEX10','CNOT4','ARID4B','CLPX','CHIC2','CXCL2','FBXO11','MTF2','CDK2','DNTTIP2','GADD45A','GOLT1B','POLR2K','NFKBIE','GABPB1','ECD','PHKG2','RAD9A','NET1','KIAA0753','EZH2','NRAS','ATP6V0B','CDK7','CCNH','SENP6','TIPARP','FOS','ARPP19','TFAP2A','KDM5B','NPC1','TP53BP2','NUSAP1']
dn_genes = ['SCCPDH','KIF20A','FZD7','USP22','PIP4K2B','CRYZ','GNB5','EIF4EBP1','PHGDH','RRAGA','SLC25A46','RPA1','HADH','DAG1','RPIA','P4HA2','MACF1','TMEM97','MPZL1','PSMG1','PLK1','SLC37A4','GLRX','CBR3','PRSS23','NUDCD3','CDC20','KIAA0528','NIPSNAP1','TRAM2','STUB1','DERA','MTHFD2','BLVRA','IARS2','LIPA','PGM1','CNDP2','BNIP3','CTSL1','CDC25B','HSPA8','EPRS','PAX8','SACM1L','HOXA5','TLE1','PYGL','TUBB6','LOXL1']
payload = {
	'up_genes': up_genes,
	'dn_genes': dn_genes,
	'direction': 'similar'
}
r = requests.post(BASE_URL + 'search', json=payload)
print r.status_code
response = r.json()
print response[0].keys()

## 3. retrieve a signature using id
uid = response[0]['id']
r = requests.get(BASE_URL + 'api', params={'id': uid})
sig = r.json()
print sig.keys()
