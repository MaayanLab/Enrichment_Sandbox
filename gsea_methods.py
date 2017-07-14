import csv
import os
import random
import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import scipy.stats as stats
import math
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, \
	RandomTreesEmbedding, AdaBoostClassifier, ExtraTreesClassifier
import json, requests
import operator, time

def Control(l_tf_genes, target, f_matrix_cols):
	'''Return the tfs in random order.'''
	return random.sample(list(f_matrix_cols), len(f_matrix_cols))

def Fisher(l_tf_genes, target, f_matrix):
	'''Return the tfs with ascending p vals as ranked by Fisher exact test with greater alternative.'''
	p = pd.Series(index=f_matrix.columns)
	for column in f_matrix:
		f_tf_genes = set(f_matrix[f_matrix[column]].index)
		a = len(f_tf_genes & l_tf_genes)
		b = len(f_tf_genes) - a
		c = len(l_tf_genes) - a
		d = 20000 - a - b - c
		o,p[column] =  stats.fisher_exact([[a,b],[c,d]], alternative='greater')
	p.sort_values(ascending=True, inplace=True)
	return list(p.index)

def FisherAdjusted(l_tf_genes, target, f_matrix, l_lib, f_lib):
	'''Like Fisher(), but weighs p vals by gene correlation within the intersection cell of the contingency table,
	Reward for high correlation. Also, weigh this reward by the degree of overlap with the ARCHS4 library.'''
	cwd = os.getcwd()
	os.chdir('..')
	os.chdir('libs')
	#Get the correlation data
	ARCHS4 = h5py.File(l_lib + '_ARCHS4_corr.h5', 'r+')
	os.chdir('..')
	os.chdir(cwd)
	if f_lib == 'ENCODE': organism_dict = {'hg19': 'human', 'mm9': mouse}
	#tfs within f_matrix will either be from human or mouse. So, store both gene lists in ARCHS4_genes_dict.
	#Note that the /index/ contains the genes, and the /values/ contain the indices. 
	#This is because we will need to access the indices from a list of genes. 
	h_genes = ARCHS4['human']['meta']['genes']
	ARCHS4_genes_dict = {'human': pd.Series(np.arange(len(h_genes)), index=h_genes[...])}
	if f_lib != 'ENCODE_TF_ChIP-seq_2015':
		m_genes = ARCHS4['mouse']['meta']['genes']
		ARCHS4_genes_dict['mouse'] = pd.Series(np.arange(len(m_genes)), index=m_genes[...])

	#For each tf, store the information collected in 'info' for use in the next iteration.
	info = pd.DataFrame(index=['a', 'p', 'o_frac', 'r', 'p_adjusted'], columns = f_matrix.columns)
	info.loc['o_frac',:] = 0
	for tf in list(f_matrix.columns):
		#Get the regular p val, just as we do in Fisher().
		f_genes = set(f_matrix[f_matrix[tf]].index)
		#'a_genes' are the genes in both the feature library tf and the label library tf.
		#In other words, 'a_genes' is the intersection cell of the 2x2 contingency table. 
		a_genes = f_genes & l_tf_genes
		a = len(a_genes)
		info.at['a',tf] = a
		b = len(f_genes) - info.at['a',tf]
		c = len(l_tf_genes) - info.at['a',tf]
		d = 20000 - info.at['a',tf] - b - c
		o, info.at['p', tf] = stats.fisher_exact([[a,b],[c,d]], alternative='greater')

		#Determine which organism this tf data came from. Doing this depends on the gmt file. 
		if f_lib == 'CREEDS': organism = tf.partition(' GSE')[0].rpartition(' ')[2]
		elif f_lib == 'ChEA_2016': organism = tf.rpartition('_')[2].lower()
		elif f_lib == 'ENCODE_TF_ChIP-seq_2015': organism = organism_dict(tf.rpartition('_')[2].lower())
		else: print('invalid lib name!')
		if organism in ['human', 'mouse']:
			#Use 'ARCHS4_genes_dict' to get the appropriate list of ARCHS4 genes
			ARCHS4_genes = ARCHS4_genes_dict[organism]
			#'overlap' will contain a list of genes both in ARCHS4 and in the intersection cell.
			overlap = {x.encode('utf-8') for x in a_genes} & set(ARCHS4_genes.index)
			l_o = len(overlap)
			if l_o > 0: 
				#'o_frac' is the proportion of genes in the intersection cell which are also in ARCHS4.
				#Limit its value to < .95 to prevent an extreme effect when used in the adjustment formula. 
				info.at['o_frac',tf] = min(l_o / a,.95)
				if l_o > 1:
					#Get the indices of the overlapping genes, and use this to index the correlation matrix.
					overlap_indices = sorted(ARCHS4_genes[overlap])
					r_vals = ARCHS4[organism]['data']['correlation'][overlap_indices][:,overlap_indices]
					#r_vals is the correlation matrix for only the overlapping tfs. 
					#Now, get the average r value for non-diagonal i.e. pairwise entries. 
					#(Each pair is duplicated across the diagonal, but this does not affect the result.)
					info.at['r',tf] = min((np.sum(r_vals) - l_o)/(l_o*l_o-l_o),.95)
		else: print('weird organism:', organism)
	ARCHS4.close()

	#for tfs which did not have at least 2 genes in ARCHS4, set their 'r' val to the median. 
	#(We should not set to zero, since this would "punish" them in comparison to tfs with very low 'r' val.)
	r_grand_median = info.loc['r',:].median()
	if pd.isnull(r_grand_median): r_grand_median = 0
	info.loc['r',:].fillna(r_grand_median, inplace=True)

	#Get the adjusted p val for each tf. Lower adjusted p vals are still considered more significant.
	for tf in list(f_matrix.columns):
		info.at['p_adjusted',tf] = math.log(info.at['p',tf]) * ((1/(1-info.at['r',tf])) ** info.at['o_frac',tf])
	ordered_tfs = list(info.columns[info.loc['p_adjusted',:].argsort()])
	return ordered_tfs

def ZAndCombined(l_tf_genes, target, f_lib, l_name, f_tfs):
	'''Uses the Enrichr API to return two lists containing the Z score and Combined score rankings.
	Note: results are not exactly the same: my ties are in different order.'''

	def get_id(l_tf_genes):
		'''Give Enrichr the list of genes in the label tf. Returns the user_list_id.'''
		ENRICHR_URL_ID = 'http://amp.pharm.mssm.edu/Enrichr/addList'
		genes_str = '\n'.join(l_tf_genes)
		description = l_name
		payload = {
			'list': (None, genes_str),
			'description': (None, l_name)
		}

		response = requests.post(ENRICHR_URL_ID, files=payload)
		if not response.ok:
		    raise Exception('Error analyzing gene list')

		data = json.loads(response.text)
		return(data['userListId'])

	#I believe Enrichr does not always return results for low-rankings feature tfs. 
	#So, set all scores beforehand to an impossibly-low score, so that they are still returned in the rankings.
	z_scores = {k:1000000 for k in f_tfs}
	combined = {k:-1000000 for k in f_tfs}
	if f_lib == 'CREEDS': libs = ['Single_Gene_Perturbations_from_GEO_up', 'Single_Gene_Perturbations_from_GEO_down']
	else: libs = [f_lib]
	for x in libs:
		#Give Enrichr the user_list_id and gmt library name to get the enrichment results.
		ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
		query_string = '?userListId=%s&backgroundType=%s'
		gene_set_lib = x
		user_list_id = get_id(l_tf_genes)
		url = ENRICHR_URL + query_string % (user_list_id, gene_set_lib)
		print(url)
		time.sleep(1) #Delay needed, otherwise Enrichr returns an error

		response = requests.get(url)
		if not response.ok:
		    raise Exception('Error fetching enrichment results')
		data = json.loads(response.text)

		#Collect the Z and Combined scores for each tf. 
		#For CREEDS, keep the best score between ...GEO_up and ...GEO_down. 
		data = data[gene_set_lib]
		for tf in data:
			if tf[1] not in z_scores: print('ERROR: ' + tf[1] + ' is not in overlap list.')
			if z_scores[tf[1]] == 1000000:
				z_scores[tf[1]] = tf[3]
				combined[tf[1]] = tf[4]
			else:
				z_scores[tf[1]] = min(z_scores[tf[1]], tf[3])
				combined[tf[1]] = max(combined[tf[1]], tf[4])

	sorted_z_scores = sorted(z_scores.items(), key=operator.itemgetter(1), reverse=False)
	sorted_combined = sorted(combined.items(), key=operator.itemgetter(1), reverse=True)
	return list(pd.DataFrame(sorted_z_scores)[0]), list(pd.DataFrame(sorted_combined)[0])

def Forest(l_tf_genes, target, train_group, features, random_state, max_features, bootstrap, class_weight, max_depth):
	clf = RandomForestClassifier(random_state = random_state, max_features=max_features, 
		bootstrap = bootstrap, class_weight = class_weight, max_depth=max_depth)
	clf.fit(train_group[features], target)
	importances = clf.feature_importances_
	rankings = pd.Series(importances, index=features)
	rankings.sort_values(ascending=False, inplace=True)
	return pd.Series(list(rankings.index))

def ForestDrop(l_tf_genes, target, train_group, features, random_state, max_features, bootstrap, class_weight, max_depth):
	f = list(features)
	rankings = pd.Series(index=features)
	x = 0
	while x < len(list(features)):
		if x<50: n_to_drop = 1
		else: n_to_drop = 300
		this_iteration_ranks = Forest(l_tf_genes, target, train_group, f, 
                                random_state, max_features, bootstrap, class_weight, max_depth)
		top_features = this_iteration_ranks[0:n_to_drop]
		for tf in top_features:
			rankings[tf] = x
			x += 1
			f.remove(tf)
	rankings.sort_values(ascending=True, inplace=True)
	return pd.Series(list(rankings.index))

def ML_wrapper(l_tf_genes, target, method, train_group, features, random_state):
	clf = method(random_state = random_state)
	clf.fit(train_group[features], target)
	importances = clf.feature_importances_
	rankings = pd.Series(importances, index=features)
	rankings.sort_values(ascending=False, inplace=True)
	return pd.Series(list(rankings.index))