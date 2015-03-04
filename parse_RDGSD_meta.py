## To parse the meta data from RDGSD folder
import os
os.chdir('D:\Zichen_Projects\microtask_GEO\RDGSD')

with open ('RDGSD_fev01_cleaned.txt') as f:
	header = next(f).split('\t')
	print header
	# for line in f:
