import csv
import os
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import enrichment_methods as m
from setup import convert_gmt
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, \
	RandomTreesEmbedding, AdaBoostClassifier, ExtraTreesClassifier
import h5py
from setup import open_csv
from math import log

def enrichment_wrapper(pair):
	'''This function is called for each lib pair, and iterates over each method and each tf. 
	pair : dict
		key 'l' contains the label library df, and key 'f' contains the feature library df. 
	'''
	l_name, f_name = pair['l'], pair['f']
	output_heading = 'from_' + l_name + '_to_' + f_name

	fisher = open_csv(output_heading+'_Fisher.csv')
	forest = open_csv(output_heading+'_RandomForest.csv')
	CombinedFF = pd.DataFrame(index=fisher.index)
	for x in fisher.columns:
		FF = [-log(max(fi, 1e-100))*(fo + 1e-3) for (fi,fo) in zip(fisher[x], forest[x])]
		CombinedFF[x] = FF
		print(CombinedFF.shape)
	CombinedFF.to_csv(output_heading+'_CombinedFF2.csv', sep='\t')
	return

if __name__ == '__main__':
	all_libs = ['CREEDS', 'ENCODE_TF_ChIP-seq_2015', 'ChEA_2016']

	os.chdir('results')

	#Iterate over each gmt pair.
	lib_pairs = [{'l':a, 'f':b} for a in all_libs for b in all_libs if a != b]
	Parallel(n_jobs=6, verbose=0)(delayed(enrichment_wrapper)(pair)for pair in lib_pairs)