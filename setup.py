import requests
import csv
import os
import h5py
import pickle
import numpy as np
import pandas as pd
from joblib import Parallel, delayed 

def open_csv(fname):
	return pd.read_csv(fname, index_col=0, sep='\t', low_memory=False, encoding='Latin-1')

def open_gvm(fname):
	print('opening', fname)
	#Read in chunks if the file is too big.
	if 'interactions' in fname:
		file_chunks = pd.read_csv(fname, keep_default_na = False, sep='\t', 
			low_memory=False, encoding='Latin-1', index_col=0, chunksize=1000)
		gvm = pd.concat(file_chunks)
		gvm = gvm.replace(to_replace='',value=False).astype(bool)
	else:
		gvm = pd.read_csv(fname, keep_default_na = False, na_values=('',), sep='\t', 
			low_memory=False, encoding='Latin-1', index_col=0)
		##Workaroud from bug which raised error if keep_default_na=False and index_col=0 are both used.
		#gvm = pd.read_csv(fname, keep_default_na = False, na_values=('',), sep='\t', low_memory=False, encoding='Latin-1')
		#lib_name = fname.partition('_gvm.csv')[0]
		#gvm.set_index(gvm[lib_name], inplace=True)
		#gvm.drop([lib_name], axis=1, inplace=True)

		#Convert blank cells to False, and ensure bool type. 
		gvm = gvm.fillna(False).astype(bool)
	return gvm

def file_exists(fname):
	'''Checks if a file exists in the directory, printing a statement if so.'''
	if os.path.isfile(fname):
		print(fname, 'has already been created.')
		return True
	else: return False

def convert_gmt(gmt_fname, output_type='gvm'):
	'''
	Converts a gmt file to either a dict or a gvm dataframe, and then returns it.
	If converting to a gvm, the resulting dataframe is also saved as a csv file in the current working directory.
	gmt_fname : str
		The gmt filename, i.e. "ChEA_2016.txt". The file must be in the current working directory. 
	output_type : str
		Is either 'gvm' or 'dict'.
	'''
	#The first segment of this list, ending at 'nan', comes from pandas.read_csv(na_values) documentation.
		#From the original list, 'NA' was removed because it is, in fact, a gene. 
	#The second segment of this list, beginning with '---', was added according to my own observations.
	MY_NA_VALS = {'', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN',
		'-NaN', '-nan', '1.#IND', '1.#QNAN', 'N/A', 'NULL', 'NaN', 'nan', 
		'---', '[NULL]'}

	def gmt_to_gvm(reader, output_fname):
		'''
		Transforms a gmt .txt file to a gvm dataframe, saves it, and returns it.
		reader : reader object
		output_fname : str
			The output filename, i.e. "ChEA_2016_gvm.csv".
		'''
		lib_name = output_fname.partition('_gvm')[0]

		gvm = pd.DataFrame(dtype=bool)
		#Each row corresponds to an annotation (e.g. a tf or drug). So, for each annotation:
		for row in reader:
			annotation = row[0]
			print(annotation)
			if annotation not in MY_NA_VALS:
				#Obtain the gene set, removing delimiting strings and not accepting null values.
				genes = [str(x).replace(',1.0', '') for x in row[2:] if str(x).replace(',1.0', '') not in MY_NA_VALS]
				#Remove repeats. To save time, comment out this below line if gene sets are known to be unique.
				genes = set(genes)
				#Create a column gene vector for this annotation...
				vec = pd.DataFrame(True, index = genes, columns = [annotation], dtype=bool)
				#...and concatenate it with the growing dataframe. 
				gvm = pd.concat([gvm,vec], axis=1)
		gvm.index.name = lib_name
		#Save file to the current working directory, with values True and '' (instead of False, to save space).
		gvm.to_csv(output_fname, sep='\t')
		#Return the dataframe in sparse form, with values True and False.
		gvm = gvm.fillna(False).to_sparse()
		return gvm

	def gmt_to_dict(reader):
		'''
		Converts a gmt .txt file to a dict, and returns it.
		reader : reader object
		'''
		d = {}
		#Each row corresponds to an annotation (e.g. a tf or drug). So, for each annotation:
		for row in reader:
			annotation = row[0]
			if annotation not in MY_NA_VALS:
				#Obtain the gene set, removing delimiting strings and not accepting null values.
				#Store the row data in the dict with the annotation as the key and the gene set as the value.
				d[annotation] = {str(x).replace(',1.0', '') for x in row[2:] if x not in MY_NA_VALS}
		return d

	#If the gvm is requested, check to see if it has already been created. If so, simply load it and return it.
	output_fname = gmt_fname.partition('.')[0] + '_gvm.csv'
	if output_type == 'gvm' and os.path.isfile(output_fname): 
		print('will use old gvm file for', gmt_fname.partition('.')[0])
		result = open_gvm(output_fname)
		#Ensure the index name is the name of the gmt.
		result.index.name = gmt_fname.partition('.')[0]
		return result

	#Else, open the gmt file.
	print('getting', gmt_fname, 'as', output_type)
	if not os.path.isfile(gmt_fname): raise ValueError(gmt_fname, 'is not in the current working directory.')	
	with open(gmt_fname, 'r') as f:
		reader = csv.reader(f, delimiter = '\t')
		#Use either gmt_to_gvm or gmt_to_dict to create and return the desired data structure.
		if output_type == 'gvm': return gmt_to_gvm(reader, output_fname)
		elif output_type == 'dict': return gmt_to_dict(reader)
		else: raise ValueError(output_type, 'must be either gvm or dict.')

def combine_gmts(gmt_fnames, output_fname, merge_type='union'):
	'''
	Creates a gvm with merged gene sets from two gmt files with corresponding annotations in the same order.
	E.g. gmt file A has down-regulated genes, and gmt file B has up-regulated genes for the same annotations. 
	gmt_fnames : list-like
		Contains the names of the two gmt files, which must be in the current working directory. 
	output_fname : str
		Name of the output file.
	merge_type : str
		'union' or 'intersection'
	'''
	def dict_to_gvm(d, output_fname):
		df = pd.DataFrame(dtype=bool)
		#Convert each key to a gene vector column in the matrix.
		for k in d:
			vec = pd.DataFrame(True, index = d[k], columns = [k], dtype=bool)
			df = pd.concat([df,vec], axis=1)
		df.index.name = output_fname.partition('_gvm.csv')[0]
		df.to_csv(output_fname, sep='\t')
		return df

	#Exit if already completed.
	if file_exists(output_fname): return
	print('creating', output_fname)

	#Otherwise, first convert each gmt to a dict.
	print('getting gmts as dicts')
	dicts = [convert_gmt(x, 'dict') for x in gmt_fnames]

	#Then, merge the dicts according to the specified method, union or intersection.
	print('merging dicts')
	combined_dict = {}
	if merge_type == 'union':
		for k in dicts[0]: combined_dict[k] = dicts[0][k] | dicts[1][k]
	elif merge_type == 'intersection': 
		for k in dicts[0]: combined_dict[k] = dicts[0][k] & dicts[1][k]

	#Finally, convert the merged dict to a gvm matrix and return it.
	#(Also save it in the current working directory.)
	print('converting merged dict to df')
	return dict_to_gvm(combined_dict, output_fname)

def combine_paired_gvm(input_fname, output_fname, merge_type='union'):
	'''
	Creates a gvm with merged gene sets from a gvm file where columns are up/down paired, 
		i.e. (col1,col2),(col3,col4),(col5,col6).
	Odd columns must be labeled "[annotation]-up" and even columns must be labeled "[annotation]-dn".
	input_fname : str
		Name of the input gvm file, which must be in the current working directory. 
	output_fname : str
		Name of the output file. 
	merge_type : str
		'union' or 'intersection'
	'''
	#Exit if already completed.
	if file_exists(output_fname): return

	old_gvm = open_gvm(input_fname)
	new_gvm = pd.DataFrame(index=old_gvm.index)
	dn_cols = list(old_gvm.columns)[::2]
	for dn_col in dn_cols:
		col = str(dn_col).rpartition('-dn')[0]
		up_col = col + '-up'
		if up_col not in old_gvm.columns: raise ValueError(up_col)
		if merge_type == 'union': new_gvm[col] = old_gvm[dn_col] | old_gvm[up_col]
		elif merge_type == 'intersection': new_gvm[col] = old_gvm[dn_col] & old_gvm[up_col]
		else: raise ValueError('invalid merge type: ' + merge_type)

	new_gvm = new_gvm.replace(to_replace=False, value='')
	new_gvm.index.name = output_fname.partition('_gvm')[0]
	new_gvm.to_csv(output_fname,sep='\t')

if __name__ == '__main__':

	os.chdir('libs')

	#Convert all transcription factor gmts to gvms.
	gmts = (
		'ChEA_2016.txt',
		'ENCODE_TF_ChIP-seq_2015.txt',
		'Single_Gene_Perturbations_from_GEO_down.txt',
		'Single_Gene_Perturbations_from_GEO_up.txt')
	Parallel(n_jobs=3, verbose=0)(delayed(convert_gmt)(gmt, 'gvm') for gmt in gmts)
	#Get the CREEDS_tfs gvm by taking the union of CREEDS' up and down gene sets.
	combine_gmts(['Single_Gene_Perturbations_from_GEO_down.txt', 'Single_Gene_Perturbations_from_GEO_up.txt'], 
		'CREEDS_TFs_gvm.csv')

	#NOTE: Drug library gvms must be downloaded from [LINK], or obtained from
	#the DrugLib_Enrichment_Comparison repo at https://github.com/MaayanLab/DrugLib_Enrichment_Comparison/ .
	#These should be manually saved into the \libs folder. 

	os.chdir('..')