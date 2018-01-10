import csv
import os
import random
import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import scipy.stats as stats
from ast import literal_eval
from math import log
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, \
	RandomTreesEmbedding, AdaBoostClassifier, ExtraTreesClassifier
import json, requests
import operator, time
from sklearn.svm import LinearSVC
from itertools import chain

#==========================================================================================
'''
Each method MUST have the FIRST argument be input_geneset, a list-like object
	containing the genes in the input library gene vector.  
Why? See perform_enrichment in get_scores.py .
'''
#==========================================================================================

def Control(input_geneset, slib_annots):
	'''
	Assigns ranks randomly: input_geneset is unused.
	slib_annots : list-like
	'''
	return random.sample(range(len(slib_annots)), len(slib_annots))

def Fisher(input_geneset, slib_gvm):
	'''
	Assigns p-values from the Fisher exact test using the greater alternative.
	Assumes the total number of genes in the space is 20,000
		for p-value comparability between enrichment results using different libraries.
	slib_gvm : pandas.DataFrame
	'''
	pvals = pd.Series(index=slib_gvm.columns)
	input_geneset = set(input_geneset)
	for annot in slib_gvm:
		annot_geneset = set(slib_gvm.index[slib_gvm[annot]])
		a = len(annot_geneset & input_geneset)
		b = len(annot_geneset) - a
		c = len(input_geneset) - a
		d = 20000 - a - b - c
		o,pvals[annot] =  stats.fisher_exact([[a,b],[c,d]], alternative='greater')
	return(tuple(pvals))

def BinomialProportions(input_geneset, slib_gvm):
	'''
	A variation of the Fisher's exact test which uses the actual number of genes in 
		the space of the search library, instead of assumming there are 20,000. 
	slib_gvm : pandas.DataFrame
	'''
	pvals = pd.Series(index=slib_gvm.columns)
	input_geneset = set(input_geneset)
	for annot in slib_gvm:
		annot_geneset = set(slib_gvm.index[slib_gvm[annot]])
		a = len(annot_geneset & set(input_geneset))
		b = len(annot_geneset)
		c = len(input_geneset) - a
		d = slib_gvm.shape[0] - b
		o,pvals[annot] =  stats.fisher_exact([[a,b],[c,d]], alternative='greater')
	return(tuple(pvals))

def ZAndCombined(input_geneset, slib_name, slib_annots):
	'''
	Uses the Enrichr API to return two lists containing the Z score and Combined score rankings.
	Note: results are not exactly the same: my ties are in a different order.
	slib_name : str
	slib_annots : list-like
	'''
	def get_id(input_geneset):
		'''Give Enrichr the input gene set to obtain a userListId.'''
		ENRICHR_URL_ID = 'http://amp.pharm.mssm.edu/Enrichr/addList'
		genes_str = '\n'.join(input_geneset)
		payload = {
			'list': (None, genes_str),
		}

		response = requests.post(ENRICHR_URL_ID, files=payload)
		if not response.ok:
		    raise Exception('Error analyzing gene list')

		data = json.loads(response.text)
		return(data['userListId'])

	#I believe Enrichr does not always return results for low-scoring annotations. 
	#So, set all scores beforehand to an impossibly-low score, so that these annotations are still assigned a score.
	z_scores = pd.Series(1000000, index=slib_annots)
	combined = pd.Series(-100000, index=slib_annots)

	if slib_name == 'CREEDS_tfs': slibs = ('Single_Gene_Perturbations_from_GEO_up', 'Single_Gene_Perturbations_from_GEO_down')
	#elif slib_name == 'CREEDS_drugs': ...
	else: slibs = (slib_name)
	for slib in slibs:
		#Give Enrichr the user_list_id and gmt library name to get the enrichment results.
		ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
		query_string = '?userListId=%s&backgroundType=%s'
		gene_set_lib = slib
		user_list_id = get_id(input_geneset)
		url = ENRICHR_URL + query_string % (user_list_id, gene_set_lib)
		time.sleep(1) #Delay needed, otherwise Enrichr returns an error.

		response = requests.get(url)
		if not response.ok:
		    raise Exception('Error fetching enrichment results')
		data = json.loads(response.text)
		results = data[gene_set_lib]

		#Collect the Z and Combined scores for each annotation. 
		#For CREEDS, keep the best score between ...GEO_up and ...GEO_down. 
		for annot_result in results:
			annot = annot_result[1]
			if annot not in slib_annots: raise Exception('ERROR: ' + annot + ' is not a search library annotation.')
			if z_scores[annot] == 1000000:
				z_scores[annot] = annot_result[3]
				combined[annot] = annot_result[4]
			else:
				z_scores[annot] = min(z_scores[annot], annot_result[3])
				combined[annot] = max(combined[annot], annot_result[4])

	return tuple(z_scores), [-score for score in tuple(combined)]

def ML_wrapper(input_geneset, method, train_group, features, random_state):
	'''This is a wrapper for sklearn.ensemble methods.'''
	target = [str(x) in input_geneset for x in train_group.index.values]
	clf = method(random_state = random_state)
	clf.fit(train_group[features], target)
	if method == LinearSVC: return [-abs(x) for x in clf.coef_[0]]
	else: return [-score for score in clf.feature_importances_]

def ML_wrapper_2(input_geneset, method, train_group, features, random_state, max_depth):
	'''This is a duplicate of ML_wrapper, in case parameters need to be added to the method.'''
	target = [str(x) in input_geneset for x in train_group.index.values]
	clf = method(random_state = random_state, max_depth=max_depth)
	clf.fit(train_group[features], target)
	if method == LinearSVC: return [-abs(x) for x in clf.coef_[0]]
	else: return [-score for score in clf.feature_importances_]

def ML_iterative(input_geneset, method, it, train_group, features, random_state):
	'''
	This is a wrapper for sklearn.ensemble methods which calls the method recursively, incrementally choosing features.
	it : int
		the number of features chosen at each iteration, within the top fifty ranks. 
	'''
	f = list(features)
	scores = pd.Series(index=features)
	x = 0
	while x < len(list(features)):
		#For the top fifty features, only choose "it" number of features at each iteration.
		if x<50: n_to_drop = it
		else: n_to_drop = 300
		this_iteration_scores = pd.Series(ML_wrapper(input_geneset, method, train_group, f, random_state), index = f)
		top_features = list(this_iteration_scores.sort_values().index)[0:n_to_drop]
		#Take the best features, then call the method again using all but the best features.
		for tf in top_features:
			scores[tf] = -1000000 + x
			x += 1
			f.remove(tf)
	return scores

def ML_fisher_cutoff(input_geneset, method, cutoff_frac, train_group, features, random_state):
	'''Uses Fisher to remove low-ranking features, then uses ML on the rest to rank them.'''
	fisher_results = pd.Series(Fisher(input_geneset, train_group), index=train_group.columns)
	n_to_keep = int(len(fisher_results) * cutoff_frac)
	new_features = fisher_results.sort_values().index[:n_to_keep]
	new_train_group = train_group[new_features]
	#Call the ML method using the remaining subset.
	new_scores_for_top_features = pd.Series(ML_wrapper(input_geneset, RandomForestClassifier, new_train_group, new_features, random_state), 
		index=new_features)
	#range(Fisher) is [0,1]; range(ML_wrapper) is [-1,0] so new_scores_for_top_features are guaranteed to be in a higher rank tier, so to speak.
	for f in new_features: fisher_results[f] = new_scores_for_top_features[f]
	return fisher_results.values

def ML_fisher_cutoff_V2(input_geneset, method, cutoff_frac, train_group, features, random_state):
	'''Uses Fisher to remove low-ranking features, then uses ML on all features to rank only the remaining subset.'''
	fisher_results = pd.Series(Fisher(input_geneset, train_group), index=train_group.columns)
	n_to_keep = int(len(fisher_results) * cutoff_frac)
	top_features = fisher_results.sort_values().index[:n_to_keep]
	#Call the ML method using all of the features...
	ML_scores_for_top_features = pd.Series(ML_wrapper(input_geneset, RandomForestClassifier, train_group, features, random_state), 
		index=features)
	#...but only revise scores for the subset which ranked highly on the Fisher test.
	for f in top_features: fisher_results[f] = ML_scores_for_top_features[f]
	return fisher_results.values

def Fisher_ML_cutoff(input_geneset, method, cutoff_frac, train_group, features, random_state):
	'''Uses ML to remove low-ranking features, then uses RandomForest to rank the rest.'''
	ML_results = pd.Series(ML_wrapper(input_geneset, RandomForestClassifier, train_group, features, random_state), index=features)
	n_to_keep = int(len(ML_results) * cutoff_frac)
	top_features = ML_results.sort_values().index[:n_to_keep]
	#Get p values for the top features
	p_vals_for_top_features= pd.Series(Fisher(input_geneset, train_group[top_features]), index=top_features)
	#Revise the scores of the top features. range(Fisher - 2) will be [-2,0], which sets it in a higher rank tier, so to speak.
	for f in top_features: ML_results[f] = p_vals_for_top_features[f] - 2
	return ML_results.values

def Impurity(input_geneset, slib_gvm, metric):
	'''
	Returns the annotations with ascending impurity.
	metric : str
		Either "Entropy" or "Gini"
	'''

	def I(p,n, metric):
		a, b = p/(p+n), n/(p+n)
		#Necessary to check if p or n is zero, because log(0) is undefined.
		if 0 in [p,n]: return 0
		elif metric == 'Entropy': return - a * log(a,2) - b * log(b,2)
		elif metric == 'Gini': return a * b
		else: raise Exception('Invalid metric.')

	#We will store the score for each input library annotation in `results`.
	results = pd.Series(index=slib_gvm.columns)
	#Get the set of genes corresponding to the input library tf.
	input_geneset = set(input_geneset)
	#Get the set of genes in all of the search library.
	slib_genes = set(slib_gvm.index)

	#Classify the search library genes as either in the input gene set (p), or not (n).
	#(This is the notation from the 1975 ID3 paper.)
	p_set = slib_genes & input_geneset
	n_set = slib_genes - p_set
	p, n = len(p_set), len(n_set)

	#For each search library annotation, calculate the info gain resulting from a split on the annotation.
	for annot in slib_gvm:
		annot_geneset = set(slib_gvm.index[slib_gvm[annot]])
		p_in_set = annot_geneset & input_geneset 
		n_in_set = annot_geneset - p_in_set

		#"in" means in the annotation gene set; "out" means not. 
		p_in, n_in = len(p_in_set), len(n_in_set)
		p_out = p - p_in
		n_out = n - n_in

		'''
		p_in is the same as "a" from Fisher.
		n_in is the same as "b" from Fisher.
		p_out is NOT the same as "c" from Fisher: it is the subset of "c" that intersects with 
			the set of genes from the search library. 
		n_out is NOT the same as "d" from Fisher. Fisher's "d" = 20000 - a - b - c.
			n_out is the genes in the search library which are not in either the input or annotation gene sets.  
		'''

		#Because the starting impurity is same for each split, it is sufficient to use the weighted average
			#for the two resulting nodes as the score.
		results[annot] = ((p_in + n_in) * I(p_in, n_in, metric) + (p_out + n_out) * I(p_out, n_out, metric)) / (p + n)
	return(list(results))

#REMOVED-DOES NOT WORK. see get_classifiers.py in the old or unused scripts folder.
# def ML_Fisher_features(input_geneset, slib_gvm, classifier):
# 	'''
# 	Applies machine learning to classify annotation gene set pairs based on the Fisher contingency table and p value.
# 	See get_classifiers in get_scores.py .
# 	'''

# 	#Build test_df.
# 	test_df = pd.DataFrame(columns=['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p'])
# 	input_geneset = set(input_geneset)
# 	for column in slib_gvm:
# 		annot_geneset = set(slib_gvm.index[slib_gvm[column]])
# 		a = len(annot_geneset & input_geneset)
# 		b = len(annot_geneset) - a
# 		c = len(input_geneset) - a
# 		d = 20000 - a - b - c
# 		o,p =  stats.fisher_exact([[a,b],[c,d]], alternative='greater')
# 		result = pd.Series([a,b,c,d,a/d,a/(a+b),a/(a+c),p],
# 							index=['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p'],
# 							name = column)
# 		test_df = test_df.append(result)

# 	#Use our classifier built in get_scores.py to predict the classifications of test_df. 
# 	return [a for (a,b) in classifier.predict_proba(test_df)]

# #REMOVED-DOES NOT WORK. see get_classifiers.py in the old or unused scripts folder.
# def ML_Fisher_features_3(input_geneset, slib_gvm, classifier, RFC, XGB, features, random_state):
# 	'''
# 	Variant of ML_Fisher_features which also considers random forest feature importance, 
# 		Gini impurity, XGBoost feature importance (though this is slow), and possibly other 
# 		engineered features. 
# 	'''
# 	test_df = pd.DataFrame(columns=['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p','rf','gini','xgb'])

# 	#Fill fisher's contingency table and p value columns.
# 	input_geneset = set(input_geneset)
# 	for column in slib_gvm:
# 		annot_geneset = set(slib_gvm.index[slib_gvm[column]])
# 		a = len(annot_geneset & input_geneset)
# 		b = len(annot_geneset) - a
# 		c = len(input_geneset) - a
# 		d = 20000 - a - b - c
# 		o,p =  stats.fisher_exact([[a,b],[c,d]], alternative='greater')
# 		result = pd.Series([a,b,c,d,a/d,a/(a+b),a/(a+c),p,None,None,None],
# 							index=['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p','rf','gini','xgb'],
# 							name = column)
# 		test_df = test_df.append(result)

# 	#Fill random forest feature importance column.
# 	target = [str(x) in input_geneset for x in slib_gvm.index.values]
# 	rf_clf = RFC(random_state = random_state)
# 	rf_clf.fit(slib_gvm[features], target)
# 	test_df['rf'] = [-score for score in rf_clf.feature_importances_]

# 	#Fill gini impurity column.
# 	test_df['gini'] = Impurity(input_geneset, slib_gvm, 'Gini')

# 	# Comment the below section out if you do not plan to use xgboost. 
# 	# (It takes a while to build the classifier.)
# 	# print('xgb')
# 	# xgb_clf = XGB(random_state = 73017)
# 	# xgb_clf.fit(slib_gvm[features], target)
# 	# test_df['xgb'] = [-score for score in xgb_clf.feature_importances_]

# 	#If xgboost is not used, remember to drop it from the columns here. 
# 	return [a for (a,b) in classifier.predict_proba(test_df.drop('xgb', axis=1))]

def pairwise_impurity(input_geneset, slib_gvm, metric):
	'''
	Calculates impurity for each possible pair of annotations.
	An annotation's score is an adjustable function of the ranks of its top and median scores.
	metric : str
		Either "Entropy" or "Gini"
	'''

	def I(p,n, metric):
		if 0 in (p,n): return 0
		a, b = p/(p+n), n/(p+n)
		if metric == 'Entropy': return - a * log(a,2) - b * log(b,2)
		elif metric == 'Gini': return a * b
		else: raise Exception('Invalid metric.')

	def split(c_and_i, c_not_i, j_set):
		'''Split a set of genes using sets i and j.'''
		cset1 = len((c_and_i) & j_set) #i and j
		cset2 = len(c_and_i) - cset1 #i not j
		cset3 = len((c_not_i) & j_set) #j not i
		cset4 = len(c_not_i) - cset3 #neither i nor j
		return (cset1, cset2, cset3, cset4)

	def final_score_function(scores, function):
		if function == '1': return np.mean(scores)
		elif function == '2': return np.mean(scores[:int(len(scores)/5)])
		elif function == '3': return np.mean(scores[:int(len(scores)/10)])
		elif function == '4': return np.mean(scores[:int(len(scores)/25)])

	#Get all possible unique pairs of annotations.
	pairs = [(i,j) for i in slib_gvm for j in slib_gvm if str(i) > str(j)]
	
	#We will store the score for each annotation pair in `pair_results`.
	pair_results = pd.Series(index=pd.MultiIndex.from_tuples(pairs, names=['tf1','tf2']))

	#We will also store the score for each individual annotation. 
	#Since there are four possible final score functions, we need four `Series`.
	individual_results1 = pd.Series(index=slib_gvm.columns)
	individual_results2 = pd.Series(index=slib_gvm.columns)
	individual_results3 = pd.Series(index=slib_gvm.columns)
	individual_results4 = pd.Series(index=slib_gvm.columns)

	#Get the set of genes corresponding to the input library annotation.
	input_geneset = set(input_geneset)

	#Get the set of genes in all of the search library.
	slib_genes = set(slib_gvm.index)

	#Classify the search library genes as either in the input gene set (p), or not (n).
	#(This is the notation from the 1975 ID3 paper.)
	p_set = slib_genes & input_geneset
	n_set = slib_genes - p_set

	#Pre-compute the sets for each annotation in slib_gvm
	annot_dict = {annot:set(slib_gvm.index[slib_gvm[annot]]) for annot in slib_gvm}

	#For each search library annotation pair, calculate the info gain resulting from the split.
	#(The split will produce four subsets whose union is the original sample, input_geneset).
	for i in set(pair_results.index.get_level_values('tf1')): #iterate over all unique annotations
		i_set = annot_dict[i]
		p_and_i, n_and_i = p_set & i_set, n_set & i_set
		p_not_i, n_not_i = p_set - p_and_i, n_set - n_and_i

		for j in pair_results[i].index:
			j_set = annot_dict[j]

			#Split the input gene set using sets i and j.
			#Obtain the sizes of the resulting sets.
			#This is like obtaining the cell counts from a 2x2 contingency table. 
			psets = split(p_and_i, p_not_i, j_set)

			#Do the same for the genes in the search library but not in the input gene set. 
			nsets = split(n_and_i, n_not_i, j_set)

			#Pair the corresponding sets from p and n. 
			sets = zip(psets, nsets)

			#Calculate the impurity score and store its value. 
			pair_results[i,j] = sum([(k[0] + k[1]) * I(k[0], k[1], metric) for k in sets])
			
	pair_results.sort_values(inplace=True)

	for annot in slib_gvm.columns:
		#Get the pairs in which the annotation appears first.
		try:aggregated_pair_scores = list(pair_results[annot])
		except KeyError: aggregated_pair_scores = []

		#Add to that the pairs in which it appears second.
		try: aggregated_pair_scores += list(pair_results[:,annot])
		except KeyError: pass

		individual_results1[annot] = final_score_function(aggregated_pair_scores, '1')
		individual_results2[annot] = final_score_function(aggregated_pair_scores, '2')
		individual_results3[annot] = final_score_function(aggregated_pair_scores, '3')
		individual_results4[annot] = final_score_function(aggregated_pair_scores, '4')

	return list(individual_results1), list(individual_results2), list(individual_results3), list(individual_results4)