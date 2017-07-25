import requests
import csv
import os
import h5py
import numpy as np
import pandas as pd
from joblib import Parallel, delayed 

def file_exists(f_name):
	'''Checks if a file exists in the directory, printing a statement if so.'''
	if os.path.isfile(f_name):
		print(f_name, 'has already been created.')
		return True
	else: return False

def convert_gmt(output_type, lib_name):
	'''
	Converts a gmt .txt file from Enrichr to either a dict or dataframe, and then returns it.
	output_type : str
		Is either 'dataframe' or 'dict'.
	lib_name : str
		Name of the gmt file. This file must be in the current working directory. 
	'''
	#The first part of this list, ending at 'nan', comes from pandas.read_csv(na_values) documentation.
		#From this list, 'NA' is removed because it is, in fact, a gene. 
	#The second part of this list, beginning with '---', was added according to my own observations.
	MY_NA_VALS = {'', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN',
		'-NaN', '-nan', '1.#IND', '1.#QNAN', 'N/A', 'NULL', 'NaN', 'nan', 
		'---', '[NULL]'}

	def to_df(reader, lib_name):
		output_fname = lib_name + '_transformed.csv'

		df = pd.DataFrame(dtype=bool)
		for row in reader:
			tf = row[0]
			print(tf)
			if tf not in MY_NA_VALS:
				#Remove delimiting string and do not accept null values
				genes = [str(x).replace(',1.0', '') for x in row[2:] if str(x).replace(',1.0', '') not in MY_NA_VALS]
				s = pd.DataFrame(True, index = genes, columns = [tf], dtype=bool)
				df = pd.concat([df,s], axis=1)
		df.index.name = lib_name
		df.to_csv(output_fname, sep='\t')
		df = df.fillna(False).to_sparse()
		return df

	def to_dict(reader, lib_name):
		d = {}
		for row in reader:
			tf = row[0]
			if tf not in MY_NA_VALS:
				d[tf] = [str(x).replace(',1.0', '') for x in row[2:] if x not in MY_NA_VALS]
		return d

	#If you want the df, check to see if it has already been created. If so, simply load it and return it.
	if output_type == 'df' and os.path.isfile(lib_name + '_transformed.csv'): 
		print('will use old df file for', lib_name)
		result = pd.read_csv(lib_name + '_transformed.csv', sep='\t', index_col=0, low_memory=False)
		result.index.name = lib_name
		return result.fillna(False).astype(bool)

	#Else, get the gmt file...
	print('getting', lib_name, 'as', output_type)
	if os.path.isfile(lib_name + '.txt'): gmt_fname = lib_name + '.txt'
	elif os.path.isfile(lib_name + '.gmt'): gmt_fname = lib_name + '.gmt'
	else: raise ValueError('No .gmt or .txt file for', lib_name)

	#And open it. Then, use either to_df or to_dict to get the data structure you want. 
	with open(gmt_fname, 'r') as f:
		reader = csv.reader(f, delimiter = '\t')
		if output_type == 'df': return to_df(reader, lib_name)
		elif output_type == 'dict': return to_dict(reader, lib_name)
		else: raise ValueError(output_type, 'must be either df or dict.')

def combine_gmts(gmts, output_fname):
	'''
	From two gmt files with the same tfs in the same order, creates a df of their gene set intersections.
	gmts : list 
		Contains the names of the two gmt files, which must be in the current working directory. 
	output_fname: str
		Name of the output file. 
	'''
	if file_exists(output_fname): return
	print('creating', output_fname)

	dicts = [convert_gmt('dict', x) for x in gmts]

	#Combine the dicts into a single dict. Note: The libraries have the same tf factors.
	combined = {}
	for k in dicts[0]: combined[k] = list(set(dicts[0][k]) | set(dicts[1][k]))

	#Convert the dict to a dataframe, one key at a time.
	df = pd.DataFrame(dtype=bool)
	for k in combined:
		s = pd.DataFrame(True, index = combined[k], columns = [k], dtype=bool)
		df = pd.concat([df,s], axis=1)
	df.to_csv('CREEDS_transformed.csv', sep='\t')
	return

def download_file(url, output_fname):
	if file_exists(output_fname): return
	r = requests.get(url, stream=True)
	with open(output_fname, 'wb') as f:
		for chunk in r.iter_content(chunk_size=1024): 
			if chunk: f.write(chunk)
	return 

def get_ARCHS4_correlation_matrices(lib):
	'''Creates .h5 files with correlation matrices between the intersection of genes from ARCHS4 and imput gmt lib files.
	libs : list
		Contains the names of the gmt files.
	'''
	new_fname = lib + '_ARCHS4_corr.h5'
	#r+ instead of w in case the operation was interrupted before, and the file was already partially created. 
	new_file = h5py.File(new_fname, 'r+')
	lib_genes = set(pd.read_csv(lib + '_transformed.csv', sep='\t', index_col=0).index)
	#Create a correlation matrix for both human genes and mouse genes. 
	for organism in ['human', 'mouse']:
		#In case the file was already partially created.
		if organism in list(new_file[...]): continue

		print(lib, organism)
		ARCHS4 = h5py.File(organism + '_matrix.h5', 'r')
		#Note the index contains the gene names. This is because we will use them to obtain the indices.
		ARCHS4_genes = pd.Series(range(len(ARCHS4['meta']['genes'])), index=ARCHS4['meta']['genes'])
		print(len(ARCHS4_genes))

		#Get the overlapping genes, and use them to get the overlapping indices.
		overlap_genes = {str(x).encode() for x in lib_genes} & set(ARCHS4_genes.index)
		print(len(overlap_genes))
		overlap_indices = ARCHS4_genes[overlap_genes].sort_values()

		#Get a sub-matrix by indexing only on the genes which were also in the gmt file
		data = pd.DataFrame(ARCHS4['data']['expression'][overlap_indices,:], index=overlap_indices.index)
		print('got data')
		data = data.transpose()
		print(data.shape)

		#Remove columns with all zero values - their Pearson correlation coefficients would be undefined. 
		data = data.loc[:, (data != 0).any(axis=0)]
		print(data.shape)
		
		#Save the genes to the .h5 file.
		genes = new_file.create_dataset(organism + '/meta/genes', data = list(data.columns))
		print('got genes')

		#Obtain and save the correlation matrix to the .h5 file.
		R = data.corr()
		print('got R', R.shape)
		corr_matrix = new_file.create_dataset(organism + '/data/correlation', data = R.values)
		print('saved R')
		
		ARCHS4.close()
	new_file.close()

if __name__ == '__main__':

	os.chdir('libs')

	#combine_gmts(['Single_Gene_Perturbations_from_GEO_down', 'Single_Gene_Perturbations_from_GEO_up'], 'CREEDS_transformed.csv')

	#This is only necessary if you plan to use the FisherAdjusted function as one of your methods.
	#for fname in ('human_matrix.h5', 'mouse_matrix.h5'):
	#	print('downloading', fname + ' . (This will take at least ten minutes.)')
	#	download_file('https://s3.amazonaws.com/mssm-seq-matrix/' + fname, fname)

	Parallel(n_jobs=2, verbose=0)(delayed(convert_gmt)('df',x) for x in ['ChEA_2016', 'ENCODE_TF_ChIP-seq_2015'])

	os.chdir('..')