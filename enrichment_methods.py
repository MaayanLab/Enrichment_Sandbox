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

#Each method MUST have the FIRST argument be l_tf_genes, the genes in the label library tf. 
#See enrichment_wrapper in get_rankings.py 

def Control(l_tf_genes, f_tfs):
	'''Return the tfs in random order.'''
	return random.sample(range(len(f_tfs)), len(f_tfs))

def Fisher(l_tf_genes, f_matrix):
	'''Return the tfs with ascending p vals as ranked by Fisher exact test with greater alternative.'''
	p = pd.Series(index=f_matrix.columns)
	for column in f_matrix:
		f_tf_genes = set(f_matrix.index[f_matrix[column]])
		a = len(f_tf_genes & set(l_tf_genes))
		b = len(f_tf_genes) - a
		c = len(l_tf_genes) - a
		d = 20000 - a - b - c
		#print(a,b,c,d) #diagnostics
		o,p[column] =  stats.fisher_exact([[a,b],[c,d]], alternative='greater')
	return(list(p))

def FisherAdjusted(l_tf_genes, f_matrix, l_lib, f_lib, ARCHS4_genes_dict):
	'''Like Fisher(), but weighs p vals by gene correlation within the intersection cell of the contingency table,
	Reward for high correlation. Also, weigh this reward by the degree of overlap with the ARCHS4 library.'''
	
	cwd = os.getcwd()
	os.chdir('..')
	os.chdir('libs')
	#Get the correlation data
	ARCHS4 = h5py.File(l_lib + '_ARCHS4_corr.h5', 'r+')
	os.chdir('..')
	os.chdir(cwd)

	#tfs within f_matrix will either be from human or mouse. So, store both gene lists in ARCHS4_genes_dict.
	#Note that the /index/ contains the genes, and the /values/ contain the indices. 
	#This is because we will need to access the indices from a list of genes. 

	if f_lib == 'ENCODE_TF_ChIP-seq_2015': organism_dict = {'hg19': 'human', 'mm9': 'mouse'}
	#For each tf, store the information collected in 'info' for use in the next iteration.
	info = pd.DataFrame(index=['a', 'p', 'o_frac', 'r', 'p_adjusted'], columns = f_matrix.columns)
	info.loc['o_frac',:] = 0
	for tf in list(f_matrix.columns):
		#Get the regular p val, just as we do in Fisher().
		f_tf_genes = set(f_matrix.index[f_matrix[tf]])
		#'a_genes' are the genes in both the feature library tf and the label library tf.
		#In other words, 'a_genes' is the intersection cell of the 2x2 contingency table. 
		a_genes = f_tf_genes & set(l_tf_genes)
		a = len(a_genes)
		info.at['a',tf] = a
		b = len(f_tf_genes) - a
		c = len(l_tf_genes) - a
		d = 20000 - a - b - c
		info.at['p', tf] = max(1e-50,stats.fisher_exact([[a,b],[c,d]], alternative='greater')[1])

		#Determine which organism this tf data came from. Doing this depends on the gmt file. 
		if f_lib == 'CREEDS': organism = tf.partition(' GSE')[0].rpartition(' ')[2]
		elif f_lib == 'ChEA_2016': 
			organism = tf.rpartition('_')[2].lower()
			if organism in ['ovary', 'hela', 'neurons', 'gbm']: organism = 'human'
		elif f_lib == 'ENCODE_TF_ChIP-seq_2015': organism = organism_dict[tf.rpartition('_')[2].lower()]
		else: print('invalid lib name!')
		if organism == 'rat': organism = 'mouse'
		if organism in ['human', 'mouse']:
			#Use 'ARCHS4_genes_dict' to get the appropriate list of ARCHS4 genes
			ARCHS4_genes = ARCHS4_genes_dict[organism]
			#'overlap' will contain a list of genes both in ARCHS4 and in the intersection cell.
			b_overlap = {str(x).encode('utf-8') for x in f_tf_genes} & set(ARCHS4_genes.index)
			c_overlap = {str(x).encode('utf-8') for x in l_tf_genes} & set(ARCHS4_genes.index)
			b_l_o = len(b_overlap)
			c_l_o = len(c_overlap)
			l_o = len(b_overlap | c_overlap)
			if l_o > 0: 
				#'o_frac' is the proportion of genes in the intersection cell which are also in ARCHS4.
				#Limit its value to < .95 to prevent an extreme effect when used in the adjustment formula. 
				info.at['o_frac',tf] = min(l_o / (b+c),.95)
				if len(b_overlap) > 0 and len(c_overlap) > 0:
					#Get the indices of the overlapping genes, and use this to index the correlation matrix.
					b_overlap_indices = sorted(ARCHS4_genes[b_overlap])
					c_overlap_indices = sorted(ARCHS4_genes[c_overlap])
					r_vals = ARCHS4[organism]['data']['correlation'][b_overlap_indices][:,c_overlap_indices]
					#r_vals is the correlation matrix for only the overlapping tfs. 
					#Now, get the average r value for non-diagonal i.e. pairwise entries. 
					#(Each pair is duplicated across the diagonal, but this does not affect the result.)
					info.at['r',tf] = min((np.sum(r_vals))/(b_l_o*c_l_o),.95)
		else: print('weird organism:', organism, tf)
	ARCHS4.close()

	#for tfs which did not have at least 2 genes in ARCHS4, set their 'r' val to the median. 
	#(We should not set to zero, since this would "punish" them in comparison to tfs with very low 'r' val.)
	r_grand_median = info.loc['r',:].median()
	if pd.isnull(r_grand_median): r_grand_median = 0
	info.loc['r',:].fillna(r_grand_median, inplace=True)

	#Get the adjusted p val for each tf. Lower adjusted p vals are still considered more significant.
	for tf in list(f_matrix.columns):
		info.at['p_adjusted1',tf] = math.log(info.at['p',tf]) * (1/(1-info.at['r',tf]) ** (5 * info.at['o_frac',tf]))
		info.at['p_adjusted2',tf] = math.log(info.at['p',tf]) * (1/(1-info.at['r',tf]) ** (10 * info.at['o_frac',tf]))
		info.at['p_adjusted3',tf] = math.log(info.at['p',tf]) - info.at['r',tf] ** (info.at['o_frac',tf])
		info.at['p_adjusted4',tf] = math.log(info.at['p',tf]) * (1 + info.at['r',tf] ** (info.at['o_frac',tf]))
		info.at['p_adjusted5',tf] = math.log(info.at['p',tf]) - (1/(1-info.at['r',tf])) ** (info.at['o_frac',tf])

	return info.loc['p_adjusted1',:], info.loc['p_adjusted2',:], info.loc['p_adjusted3',:], info.loc['p_adjusted4',:], info.loc['p_adjusted5',:]

def ZAndCombined(l_tf_genes, f_lib, f_tfs):
	'''Uses the Enrichr API to return two lists containing the Z score and Combined score rankings.
	Note: results are not exactly the same: my ties are in different order.'''
	def get_id(l_tf_genes):
		'''Give Enrichr the list of genes in the label tf. Returns the user_list_id.'''
		ENRICHR_URL_ID = 'http://amp.pharm.mssm.edu/Enrichr/addList'
		genes_str = '\n'.join(l_tf_genes)
		payload = {
			'list': (None, genes_str),
		}

		response = requests.post(ENRICHR_URL_ID, files=payload)
		if not response.ok:
		    raise Exception('Error analyzing gene list')

		data = json.loads(response.text)
		return(data['userListId'])

	#I believe Enrichr does not always return results for low-rankings feature tfs. 
	#So, set all scores beforehand to an impossibly-low score, so that they are still returned in the rankings.
	z_scores = pd.Series(1000000, index=f_tfs)
	combined = pd.Series(-100000, index=f_tfs)

	if f_lib == 'CREEDS': libs = ['Single_Gene_Perturbations_from_GEO_up', 'Single_Gene_Perturbations_from_GEO_down']
	else: libs = [f_lib]
	for x in libs:
		#Give Enrichr the user_list_id and gmt library name to get the enrichment results.
		ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
		query_string = '?userListId=%s&backgroundType=%s'
		gene_set_lib = x
		user_list_id = get_id(l_tf_genes)
		url = ENRICHR_URL + query_string % (user_list_id, gene_set_lib)
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

	return list(z_scores), [-x for x in list(combined)]

def ML_wrapper(l_tf_genes, method, train_group, features, random_state):
	#This is a wrapper for sklearn.ensemble methods.
	target = [str(x) in l_tf_genes for x in train_group.index.values]
	clf = method(random_state = random_state)
	clf.fit(train_group[features], target)
	return [-x for x in clf.feature_importances_]

def ML_iterative(l_tf_genes, method, train_group, features, random_state):
	#This is a wrapper for sklearn.ensemble methods, which chooses features recursively.
	f = list(features)
	rankings = pd.Series(index=features)
	x = 0
	while x < len(list(features)):
		if x<50: n_to_drop = 1
		else: n_to_drop = 300
		this_iteration_ranks = pd.Series(method(l_tf_genes, train_group, f, random_state), index = f)
		top_features = list(this_iteration_ranks.sort_values().index)[0:n_to_drop]
		#Take the best features, then call the method again using all but the best features.
		for tf in top_features:
			rankings[tf] = -1000000 + x
			x += 1
			f.remove(tf)
	return rankings
