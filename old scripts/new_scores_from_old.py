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

#=====================================================================================================================================
'''
The purpose of this script is to generate new score files from old score files. 
This is more efficient than running a whole new method, since the old scores have already been calculated!
In this case, a new score is being generated from the Fisher p value and RandomForestClassifier feature importance. 
i.e. New Enrichment Score = function(Fisher Score, RandomForestClassifier Score)
This script can be adapted to change the forumla used to combine the old scores.
It can also be adapted to consider different score files other than Fisher and RandomForestClassifier - even more than two at a time. 
'''
#=====================================================================================================================================

def enrichment_wrapper(pair):
	'''This function is called for each lib pair, and iterates over each tf. 
	pair : dict
		key 'l' contains the label library df, and key 'f' contains the feature library df. 
	'''
	l_name, f_name = pair['l'], pair['f']
	output_heading = 'from_' + l_name + '_to_' + f_name

	#Open the old score files.
	###This is where you can change the old scores that are being combined.###
	fisher = open_csv(output_heading+'_Fisher.csv')
	forest = open_csv(output_heading+'_RandomForest.csv')

	#Create DataFrame to store new scores.
	new_scores = pd.DataFrame(index=fisher.index)

	#Iterate over each label library tf.
	for x in fisher.columns:

		#Get the new scores, for this label library tf (the column), for each feature library tf (the rows).
		###This is where you can change the forumla used to combine the old scores.###
		new_score = [fo if fi < .9 else fo + 1 for (fi,fo) in zip(fisher[x], forest[x])]

		#Save the result to the DataFrame. 
		new_scores[x] = new_score

	new_scores.to_csv(output_heading+'_CombinedFF4.csv', sep='\t')
	return

if __name__ == '__main__':
	all_libs = ['CREEDS', 'ENCODE_TF_ChIP-seq_2015', 'ChEA_2016']

	os.chdir('results')

	#Iterate over each gmt pair.
	lib_pairs = [{'l':a, 'f':b} for a in all_libs for b in all_libs if a != b]
	Parallel(n_jobs=6, verbose=0)(delayed(enrichment_wrapper)(pair)for pair in lib_pairs)