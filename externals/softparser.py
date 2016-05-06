"""This module contains functions for parsing SOFT files.

__authors__ = "Gregory Gundersen, Andrew Rouillard, Axel Feldmann, Kevin Hu, Zichen Wang"
__credits__ = "Yan Kou, Avi Ma'ayan"
__contact__ = "gregory.gundersen@mssm.edu"
"""

from StringIO import StringIO
import numpy as np
import pandas as pd

from files import SOFTFile, ANNOTFile
from log import pprint


def parse(filename, A_cols, B_cols, annot_filename=None):
	"""Parses SOFT files, and convert probe IDs to gene sybmols if a 
	annot_filename is provided.
	"""	
	all_cols = A_cols + B_cols
	if 'GDS' in filename:
		BOF = '!dataset_table_begin'
		EOF = '!dataset_table_end'
	else:
		BOF = '!series_matrix_table_begin'
		EOF = '!series_matrix_table_end'

	# Collect lines between BOF and EOF
	lines = [] 
	with open(filename, 'r') as soft_in:
		# Skip comments.
		discard = next(soft_in)
		while discard.rstrip() != BOF:
			discard = next(soft_in)

		for line in soft_in:
			if line.rstrip() != EOF:
				lines.append(line)
	
	# Read those lines using pandas
	expr_df = pd.read_csv(StringIO(''.join(lines)), sep='\t', na_values='null')

	# Get the column name for probes
	probe_column = expr_df.columns[0]
	# Strip any potential quotation marks.
	# expr_df[probe_column] = expr_df[probe_column].apply(lambda x: x.replace('"', '').replace('\'', ''))
	
	# Subset columns
	expr_df = expr_df.set_index(probe_column)[all_cols]

	if annot_filename is not None:
		# Parse annot file
		probe2gene = parse_annot(annot_filename)
		# Make the index of expr_df to be string
		expr_df.index = expr_df.index.map(str)
		# Inner join with expr_df to get gene symbols
		expr_df = expr_df.merge(probe2gene, left_index=True, right_index=True, how='inner')
		# Set gene as index and remove probe column
		expr_df = expr_df.reset_index().set_index('gene')[all_cols]
		# Average genes 
		expr_df = expr_df.groupby(expr_df.index).aggregate(np.mean)

	# Drop genes with missing values
	expr_df = expr_df.dropna()
	return expr_df


def parse_annot(filename):
	"""Parse .annot file, return a dict of dict {platform: {probe_id: gene_symbol, ...}}
	"""
	pprint('Parsing ANNOT file.')

	probes = []
	genes = []
	BOF = '!platform_table_begin'
	EOF = '!platform_table_end'

	with open(filename, 'r') as annot_in:
		# Skip comments.
		discard = next(annot_in)
		while discard.rstrip() != BOF:
			discard = next(annot_in)

		# Read header
		header = next(annot_in).rstrip('\r\n').split('\t')

		# Find the columns indices.
		probe_index = header.index('ID')
		gene_index = header.index('Gene symbol')

		for line in annot_in:
			split_line = line.rstrip('\r\n').split('\t')
			if split_line[0] == EOF or split_line[1] == '--Control':
				continue

			probe  = split_line[probe_index]
			probe = probe.replace('"', '').replace('\'', '')
			gene = split_line[gene_index]

			# Skip probes with no corresponding genes
			if gene == '':
				continue
			# Three forward slashes, \\\, denotes multiple genes.
			if '\\\\\\' in probe:
				continue

			probes.append(probe)
			genes.append(gene)

	probe2gene = pd.DataFrame({'gene': genes, 'probe': probes}).set_index('probe')
	return probe2gene


def platform_supported(platform):
	if platform not in PROBE2GENE:
		return False
	return True


def _probe2gene(platform, probe):
	"""Converts probe IDs to gene symbols. Does not check if the platform is
	supported.
	"""

	# Strip any potential quotation marks.
	probe = probe.replace('"', '').replace('\'', '')
	try:
		if probe in PROBE2GENE[platform]:
			return PROBE2GENE[platform][probe]
	# This should never occur, given that we check if the platform is in the
	# dictionary. But just in case.
	except AttributeError:
		return None
	return None

	