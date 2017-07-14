import csv
import os
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import gsea_methods as m

def gmt_to_df(lib_name):
	'''Input: library file from Enrichr. Output: transformed library file with genes as indices and tfs as columns.
	Membership denoted by True/False.'''

	#check beforehand if the transformed file has already been created
	df_fname = lib_name + '_transformed.csv'
	if os.path.isfile(df_fname):
		print('using old df file for ' + lib_name)
		df = pd.read_csv(df_fname, index_col=0, sep='\t')
	else:
		print('transforming', lib_name)
		with open(lib_name + '.txt', 'r') as f:
			df = pd.DataFrame(False, index = [''], columns = [''], dtype=bool)
			reader = csv.reader(f, delimiter='\t')
			#reads and appends each row to the dataframe
			for row in reader:
				row = [x.replace(',1.0', '') for x in row]
				#print('transforming', row[0])
				s = pd.DataFrame(True, index = row[2:], columns = [row[0]], dtype=bool)
				df = pd.concat([df,s], axis=1)
		#try many methods of removing null values
		df = df[pd.notnull(df.index)].fillna(False)
		df = df.loc[pd.notnull(df.index)]
		df.drop('', inplace=True)
		df.drop('', axis=1, inplace=True)
		df.to_csv(df_fname, sep='\t')
		df = df.to_sparse()
	df.name = lib_name
	return df

def clean(tf):
	'''Input: a list or set of transcription factor names. Output: a set of the 'cleaned', i.e. simplified tf names.'''
	return tf.partition('_')[0].partition(' ')[0].upper()

def get_overlaps(label_tfs, feature_tfs):
	'''Input: list with feature_lib and label_lib names; list or set of feature library tf names; list or set of label library tf names.
	Output: list of tfs in label library if the clean version of that tf is also found in the feature library'''
	cleaned_overlaps = {clean(x) for x in feature_tfs} & {clean(x) for x in label_tfs}
	l_in_overlaps = [x for x in label_tfs if clean(x) in cleaned_overlaps]
	return l_in_overlaps

def get_methods(l,f, l_name, f_name):
	df = pd.DataFrame(index=['method', 'params'])
	df['Control'] = [m.Control, ([f.columns.values])]
	df['Fisher'] = [m.Fisher, ([f])]
	df['FisherAdjusted'] = [m.FisherAdjusted, (f, l_name, f_name)]
	df['ZAndCombined'] = [m.ZAndCombined, (f_name, l_name, f.columns.values)]

	train_group, features,  = f, f.columns.values, 
	random_state = 70317 #int(time.strftime('%m_%d_%y'))

	Tree = {m.Forest, m.ForestDrop}
	for x in Tree:
		n = 1
		for mf in ['auto', 'log2']:
			for bs in [True, False]:
				for cw in [None, 'balanced']:
					if x == m.ForestDrop: maxdepths = [7,5]
					else: maxdepths = [None, 7]
					for md in maxdepths:
						xname = x.__name__ + str(n)
						df[xname] = [x, (train_group, features, random_state, mf, bs, cw, md)]
						n += 1

	ML = {m.ExtraTreesClassifier, m.RandomTreesEmbedding, m.AdaBoostClassifier, 
		m.GradientBoostingClassifier}
	for x in ML:
		xname = x.__name__.partition('Classifier')[0]
		df[xname] = [m.ML_wrapper, (x, train_group, features, random_state)]

	return df

def enrichment_wrapper(l,f_index, output_dir, overlaps, m_name, mp):
	if m_name == 'ZAndCombined': 
		output_fnames = (output_dir + '_Z.csv', output_dir + '_Combined.csv')
		dfs = (pd.DataFrame(), pd.DataFrame())
	else: 
		output_fnames = (output_dir + '_' + m_name + '.csv',)
		dfs = (pd.DataFrame(),)
	if os.path.isfile(output_fnames[0]): print('ranking file already created for', m_name)
	else:
		for l_tf in list(overlaps):
			#print(m_name, l_tf)
			l_tf_genes = set(l[l[l_tf]].index)
			target = [str(x) in l_tf_genes for x in f_index]
			result = mp['method'](l_tf_genes, target, *mp['params'])
			for x in range(len(dfs)): dfs[x][l_tf] = result[x]
		for x in range(len(dfs)): dfs[x].to_csv(output_fnames[x], sep='\t')
		print('done ' + m_name)
	return

if __name__ == '__main__':
	all_libs = ['CREEDS', 'ChEA_2016'] #'ENCODE_TF_ChIP-seq_2015'

	os.chdir('libs')
	all_dfs = {x:gmt_to_df(x) for x in all_libs}
	os.chdir('..')

	lib_df_pairs = [{'l':all_dfs[a], 'f':all_dfs[b]} for a in all_dfs for b in all_dfs if a != b]
	for pair in lib_df_pairs:
		output_dir = 'from_' + pair['l'].name + '_to_' + pair['f'].name
		print('Beginning enrichment analysis. Storing results in ' + output_dir)
		if not os.path.exists(output_dir): os.makedirs(output_dir)
		os.chdir(output_dir)
		overlaps = get_overlaps(pair['l'].columns.values, pair['f'].columns.values)
		method_params = get_methods(pair['l'], pair['f'], pair['l'].name, pair['f'].name)
		#ARCHS4
		Parallel(n_jobs=6, verbose=0)(delayed(enrichment_wrapper)(pair['l'], pair['f'].index, \
			output_dir, overlaps, m_name, method_params[m_name]) for m_name in method_params)
		os.chdir('..')