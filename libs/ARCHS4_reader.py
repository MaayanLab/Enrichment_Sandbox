import h5py
import csv
import os
import numpy as np
import pandas as pd

def gmt_to_df(lib_name):
	'''Input: library file from Enrichr. Output: transformed library file with genes as indices and tfs as columns; membership denoted by True/False.
	Checks beforehand if the transformed file has already been created.'''
	df_fname = lib_name + '_transformed.csv'
	if os.path.isfile(df_fname):
		print('using old df file for ' + lib_name)
		df = pd.read_csv(df_fname, index_col=0, sep='\t')
	else:
		print('transforming', lib_name)
		with open(lib_name + '.txt', 'r') as f:
			df = pd.DataFrame(False, index = [''], columns = [''], dtype=bool)
			reader = csv.reader(f, delimiter='\t')
			for row in reader:
				row = [x.replace(',1.0', '') for x in row]
				#print('transforming', row[0])
				s = pd.DataFrame(True, index = row[2:], columns = [row[0]], dtype=bool)
				df = pd.concat([df,s], axis=1)
		df = df[pd.notnull(df.index)].fillna(False)
		df = df.loc[pd.notnull(df.index)]
		df.drop('', inplace=True)
		df.drop('', axis=1, inplace=True)
		df.to_csv(df_fname, sep='\t')
		df = df.to_sparse()
	return df

if __name__ == '__main__':
	libs = ['ENCODE_TF_ChIP-seq_2015', 'ChEA_2016', 'CREEDS']
	for lib in libs:
		new_fname = lib + '_ARCHS4_corr.h5'
		#if not os.path.exists(new_fname):
		new_file = h5py.File(new_fname, 'r+') #w
		lib_genes = set(pd.read_csv(lib + '_transformed.csv', sep='\t', index_col=0).index)
		for organism in ['human', 'mouse']:
			print(lib, organism)
			ARCHS4 = h5py.File(organism + '_matrix.h5', 'r')

			ARCHS4_genes = pd.Series(range(len(ARCHS4['meta']['genes'])), index=ARCHS4['meta']['genes'])
			print(len(ARCHS4_genes))

			overlap_genes = {str(x).encode() for x in lib_genes} & set(ARCHS4_genes.index)
			print(len(overlap_genes))

			overlap_indices = ARCHS4_genes[overlap_genes].sort_values()
			print('sorted')
			data = pd.DataFrame(ARCHS4['data']['expression'][overlap_indices,:], index=overlap_indices.index)
			print('got data')
			data = data.transpose()
			print(data.shape)

			data = data.loc[:, (data != 0).any(axis=0)]
			print(data.shape)
			
			genes = new_file.create_dataset(organism + '/meta/genes', data = list(data.columns))
			print('got genes')
			R = data.corr()
			print('got R', R.shape)
			corr_matrix = new_file.create_dataset(organism + '/data/correlation', data = R.values)
			print('saved R')
			
			ARCHS4.close()
		new_file.close()