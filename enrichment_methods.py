import csv
import os
import random
import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import scipy.stats as stats
from math import log
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, \
	RandomTreesEmbedding, AdaBoostClassifier, ExtraTreesClassifier
import json, requests
import operator, time
from sklearn.svm import LinearSVC

#==========================================================================================
'''
#Each method MUST have the FIRST argument be l_tf_genes, the genes in the label library tf. 
#See enrichment_wrapper in get_rankings.py 
'''
#==========================================================================================

def Control(l_tf_genes, f_tfs):
	'''Assigns ranks randomly.'''
	return random.sample(range(len(f_tfs)), len(f_tfs))

def Fisher(l_tf_genes, f_matrix):
	'''Assigns p vals from Fisher exact test using the greater alternative.'''
	p = pd.Series(index=f_matrix.columns)
	l_tf_genes = set(l_tf_genes)
	for column in f_matrix:
		f_tf_genes = set(f_matrix.index[f_matrix[column]])
		a = len(f_tf_genes & l_tf_genes)
		b = len(f_tf_genes) - a
		c = len(l_tf_genes) - a
		d = 20000 - a - b - c
		o,p[column] =  stats.fisher_exact([[a,b],[c,d]], alternative='greater')
	return(list(p))

def BinomialProportions(l_tf_genes, f_matrix):
	'''A variation of the Fisher's exact test.'''
	p = pd.Series(index=f_matrix.columns)
	for column in f_matrix:
		f_tf_genes = set(f_matrix.index[f_matrix[column]])
		a = len(f_tf_genes & set(l_tf_genes))
		b = len(f_tf_genes)
		c = len(l_tf_genes) - a
		d = f_matrix.shape[0] - b
		o,p[column] =  stats.fisher_exact([[a,b],[c,d]], alternative='greater')
	return(list(p))

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
		time.sleep(1) #Delay needed, otherwise Enrichr returns an error.

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
	'''This is a wrapper for sklearn.ensemble methods.'''
	target = [str(x) in l_tf_genes for x in train_group.index.values]
	clf = method(random_state = random_state)
	clf.fit(train_group[features], target)
	if method == LinearSVC: return [-abs(x) for x in clf.coef_[0]]
	else: return [-x for x in clf.feature_importances_]

def ML_wrapper_2(l_tf_genes, method, train_group, features, random_state, max_depth):
	'''This is a duplicate of ML_wrapper, in case parameters need to be added to the method.'''
	target = [str(x) in l_tf_genes for x in train_group.index.values]
	clf = method(random_state = random_state, max_depth=max_depth)
	clf.fit(train_group[features], target)
	if method == LinearSVC: return [-abs(x) for x in clf.coef_[0]]
	else: return [-x for x in clf.feature_importances_]

def ML_iterative(l_tf_genes, method, it, train_group, features, random_state):
	'''
	This is a wrapper for sklearn.ensemble methods which calls the method recursively, incrementally choosing features.
	it : int
		the number of features chosen at each iteration, within the top fifty ranks. 
	'''
	f = list(features)
	rankings = pd.Series(index=features)
	x = 0
	while x < len(list(features)):
		#For the top fifty features, only choose "it" number of features at each iteration.
		if x<50: n_to_drop = it
		else: n_to_drop = 300
		this_iteration_ranks = pd.Series(ML_wrapper(l_tf_genes, method, train_group, f, random_state), index = f)
		top_features = list(this_iteration_ranks.sort_values().index)[0:n_to_drop]
		#Take the best features, then call the method again using all but the best features.
		for tf in top_features:
			rankings[tf] = -1000000 + x
			x += 1
			f.remove(tf)
	return rankings

def ML_fisher_cutoff(l_tf_genes, method, cutoff_frac, train_group, features, random_state):
	'''Uses Fisher to remove low-ranking features, then uses ML on the rest to rank them.'''
	fisher_results = pd.Series(Fisher(l_tf_genes, train_group), index=train_group.columns)
	n_to_keep = int(len(fisher_results) * cutoff_frac)
	new_features = fisher_results.sort_values().index[:n_to_keep]
	new_train_group = train_group[new_features]
	#Call the ML method using the remaining subset.
	new_scores_for_top_features = pd.Series(ML_wrapper(l_tf_genes, RandomForestClassifier, new_train_group, new_features, random_state), 
		index=new_features)
	#range(Fisher) is [0,1]; range(ML_wrapper) is [-1,0] so new_scores_for_top_features are guaranteed to be in a higher rank tier, so to speak.
	for f in new_features: fisher_results[f] = new_scores_for_top_features[f]
	return fisher_results.values

def ML_fisher_cutoff_V2(l_tf_genes, method, cutoff_frac, train_group, features, random_state):
	'''Uses Fisher to remove low-ranking features, then uses ML on all features to rank only the remaining subset.'''
	fisher_results = pd.Series(Fisher(l_tf_genes, train_group), index=train_group.columns)
	n_to_keep = int(len(fisher_results) * cutoff_frac)
	top_features = fisher_results.sort_values().index[:n_to_keep]
	#Call the ML method using all of the features...
	ML_scores_for_top_features = pd.Series(ML_wrapper(l_tf_genes, RandomForestClassifier, train_group, features, random_state), 
		index=features)
	#...but only revise scores for the subset which ranked highly on the Fisher test.
	for f in top_features: fisher_results[f] = ML_scores_for_top_features[f]
	return fisher_results.values

def Fisher_ML_cutoff(l_tf_genes, method, cutoff_frac, train_group, features, random_state):
	'''Uses ML to remove low-ranking features, then uses Forest to rank the rest.'''
	ML_results = pd.Series(ML_wrapper(l_tf_genes, RandomForestClassifier, train_group, features, random_state), index=features)
	n_to_keep = int(len(ML_results) * cutoff_frac)
	top_features = ML_results.sort_values().index[:n_to_keep]
	#Get p values for the top features
	p_vals_for_top_features= pd.Series(Fisher(l_tf_genes, train_group[top_features]), index=top_features)
	#Revise the scores of the top features. range(Fisher - 2) will be [-2,0], which sets it in a higher rank tier, so to speak.
	for f in top_features: ML_results[f] = p_vals_for_top_features[f] - 2
	return ML_results.values

def Impurity(l_tf_genes, f_matrix, metric):
	'''Returns the tfs with ascending impurity.'''

	def I(p,n, metric):
		a, b = p/(p+n), n/(p+n)
		#Necessary to check if p or n is zero, because log(0) is undefined.
		if 0 in [p,n]: return 0
		elif metric == 'Entropy': return - a * log(a,2) - b * log(b,2)
		elif metric == 'Gini': return a * b
		else: raise Exception('Invalid metric.')

	#Store the score for each label library tf in results.
	results = pd.Series(index=f_matrix.columns)
	#Get the set of genes corresponding to the label library tf.
	l_tf_genes = set(l_tf_genes)
	#Get the set of genes in all of the feature library.
	f_lib_genes = set(f_matrix.index)

	#Classify the feature library genes as either in the label library tf set (p), or not (n).
	#(This is the notation from the 1975 ID3 paper.)
	p_set = f_lib_genes & l_tf_genes
	n_set = f_lib_genes - p_set
	p, n = len(p_set), len(n_set)

	#For each feature library tf, calculate the info gain resulting from a split on the tf.
	for column in f_matrix:
		f_tf_genes = set(f_matrix.index[f_matrix[column]])
		p_in_set = f_tf_genes & l_tf_genes 
		n_in_set = f_tf_genes - p_in_set

		#"in" means in the feature library tf set; "out" means not. 
		p_in, n_in = len(p_in_set), len(n_in_set)
		p_out = p - p_in
		n_out = n - n_in

		'''
		p_in is the same as "a" from Fisher.
		n_in is the same as "b" from Fisher.
		p_out is NOT the same as "c" from Fisher: it is the subset of "c" that intersects with 
			the set of genes from the feature library. 
		n_out is NOT the same as "d" from Fisher. Fisher's "d" = 20000 - a - b - c.
			n_out is the genes in the feature library which are not in either the label or feature tf gene sets.  
		'''

		#Because the starting impurity is same for each split, it is sufficient to use the weighted average
			#for the two resulting nodes as the score.
		results[column] = ((p_in + n_in) * I(p_in, n_in, metric) + (p_out + n_out) * I(p_out, n_out, metric)) / (p + n)
	return(list(results))

def ML_Fisher_features(l_tf_genes, f_matrix, classifier):
	'''
	Applies machine learning to classify tf gene set pairs based on Fisher contingency table and p value.
	See get_classifiers in get_scores.py .
	'''

	#Build test_df.
	test_df = pd.DataFrame(columns=['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p'])
	l_tf_genes = set(l_tf_genes)
	for column in f_matrix:
		f_tf_genes = set(f_matrix.index[f_matrix[column]])
		a = len(f_tf_genes & l_tf_genes)
		b = len(f_tf_genes) - a
		c = len(l_tf_genes) - a
		d = 20000 - a - b - c
		o,p =  stats.fisher_exact([[a,b],[c,d]], alternative='greater')
		result = pd.Series([a,b,c,d,a/d,a/(a+b),a/(a+c),p],
							index=['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p'],
							name = column)
		test_df = test_df.append(result)

	#Use our classifier built in get_scores.py to predict the classifications of test_df. 
	return [a for (a,b) in classifier.predict_proba(test_df)]

def ML_Fisher_features_3(l_tf_genes, f_matrix, classifier, RFC, XGB, features, random_state):
	'''
	Variant of ML_Fisher_features which also considers random forest feature importance, 
		Gini impurity, XGBoost feature importance (though this is slow), and possibly other 
		engineered features. 
	'''
	test_df = pd.DataFrame(columns=['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p','rf','gini','xgb'])

	#Fill fisher's contingency table and p value columns.
	l_tf_genes = set(l_tf_genes)
	for column in f_matrix:
		f_tf_genes = set(f_matrix.index[f_matrix[column]])
		a = len(f_tf_genes & l_tf_genes)
		b = len(f_tf_genes) - a
		c = len(l_tf_genes) - a
		d = 20000 - a - b - c
		o,p =  stats.fisher_exact([[a,b],[c,d]], alternative='greater')
		result = pd.Series([a,b,c,d,a/d,a/(a+b),a/(a+c),p,None,None,None],
							index=['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p','rf','gini','xgb'],
							name = column)
		test_df = test_df.append(result)

	#Fill random forest feature importance column.
	target = [str(x) in l_tf_genes for x in f_matrix.index.values]
	rf_clf = RFC(random_state = random_state)
	rf_clf.fit(f_matrix[features], target)
	test_df['rf'] = [-x for x in rf_clf.feature_importances_]

	#Fill gini impurity column.
	test_df['gini'] = Impurity(l_tf_genes, f_matrix, 'Gini')

	# Comment the below section out if you do not plan to use xgboost. 
	# (It takes a while to build the classifier.)
	# print('xgb')
	# xgb_clf = XGB(random_state = 73017)
	# xgb_clf.fit(f_matrix[features], target)
	# test_df['xgb'] = [-x for x in xgb_clf.feature_importances_]

	#If xgboost is not used, remember to drop it from the columns here. 
	return [a for (a,b) in classifier.predict_proba(test_df.drop('xgb', axis=1))]

def pairwise_impurity(l_tf_genes, f_matrix, metric):
	'''
	Calculates impurity for each possible pair of features.
	A feature's score is an adjustable function of the ranks of its top and median scores.
	'''

	def I(p,n, metric):
		if 0 in (p,n): return 0
		a, b = p/(p+n), n/(p+n)
		if metric == 'Entropy': return - a * log(a,2) - b * log(b,2)
		elif metric == 'Gini': return a * b
		else: raise Exception('Invalid metric.')

	def split(c_and_i, c_not_i, j_set):
		'''Split a set of genes using sets i and j.'''
		cset1 = (c_and_i) & j_set #i and j
		cset2 = (c_and_i) - cset1 #i not j
		cset3 = (c_not_i) & j_set #j not i
		cset4 = (c_not_i) - cset3 #neither i nor j
		return map(len, (cset1, cset2, cset3, cset4))

	def final_score_function(scores, function):
		if function == '1': return scores[0]
		elif function == '2': return np.mean(scores)
		elif function == '3': return np.mean(scores[:int(len(scores)/2)])
		elif function == '4': return np.mean([scores[0],np.mean(scores)])

	#Get all possible unique pairs of transcription factors.
	pairs = [(i,j) for i in f_matrix for j in f_matrix if str(i) > str(j)]
	
	#Store the score for each pair in ```pair_results```.
	pair_results = pd.Series(index=pd.MultiIndex.from_tuples(pairs, names=['tf1','tf2']))

	#Store the score for each transcription factor. 
	#Since there are four possible final score functions, we need four ```Series```.
	individual_results1 = pd.Series(index=f_matrix.index)
	individual_results2 = pd.Series(index=f_matrix.index)
	individual_results3 = pd.Series(index=f_matrix.index)
	individual_results4 = pd.Series(index=f_matrix.index)

	#Get the set of genes corresponding to the label library tf.
	l_tf_genes = set(l_tf_genes)

	#Get the set of genes in all of the feature library.
	f_lib_genes = set(f_matrix.index)

	#Classify the feature library genes as either in the label library tf set (p), or not (n).
	#(This is the notation from the 1975 ID3 paper.)
	p_set = f_lib_genes & l_tf_genes
	n_set = f_lib_genes - p_set

	#For each feature library tf pair, calculate the info gain resulting from the split.
	#(The split will produce four subsets whose union is the original sample, l_tf_genes).
	for i in set(pair_results.index.get_level_values('tf1')):
		print(i)
		i_set = set(f_matrix.index[f_matrix[i]])
		p_and_i, n_and_i = p_set & i_set, n_set & i_set
		p_not_i, n_not_i = p_set - p_and_i, n_set - n_and_i

		for j in pair_results[i].index:
			j_set = set(f_matrix.index[f_matrix[j]])

			#Split the label library tf set using sets i and j.
			#Obtain the sizes of the resulting sets.
			#This is like obtaining the cell counts from a 2x2 contingency table. 
			psets = split(p_and_i, p_not_i, j_set)

			#Do the same for the genes in the label library not in the tf set. 
			nsets = split(n_and_i, n_not_i, j_set)

			#Pair the corresponding sets from p and n. 
			sets = zip(psets, nsets)

			#Calculate the impurity score and store its value. 
			pair_results[i,j] = sum([(k[0] + k[1]) * I(k[0], k[1], metric) for k in sets])/(f_matrix.shape[0])

	pair_results.sort_values(inplace=True)

	for tf in f_matrix.columns:
		print(tf)
		aggregated_pair_scores = list(pair_results[tf]) + list(pair_results[:,tf])
		individual_results1[tf] = final_score_function(aggregated_pair_scores, '1')
		individual_results2[tf] = final_score_function(aggregated_pair_scores, '2')
		individual_results3[tf] = final_score_function(aggregated_pair_scores, '3')
		individual_results4[tf] = final_score_function(aggregated_pair_scores, '4')

	return list(individual_results1), list(individual_results2), list(individual_results3), list(individual_results4)