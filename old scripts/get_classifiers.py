import csv
import os
import numpy as np
import pandas as pd
from setup import convert_gmt
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, \
	RandomTreesEmbedding, AdaBoostClassifier, ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import LinearSVC
from xgboost import XGBClassifier
from setup import open_csv
import scipy.stats as stats

def get_classifiers(old_l_name, old_f_name):
	'''
	Returns classifiers which predict whether a tf gene set pair is a match or not.
	Uses ChEA in place of ENCODE, and ENCODE in place of ChEA as the training set. 
	For example, if the classifier will be used to predict enrichment from ENCODE to CREEDS,
		this function will train the classifier using enrichment from ChEA to CREEDS.
	Currently, this does not work for drug-to-gene libraries.
	old_l/f_name : str
		Names of gmt files for which the classifier will be tested on.
	'''

	#Swap ChEA for ENCODE, and ENCODE for ChEA. 
	os.chdir('..')
	os.chdir('libs')
	l_name = old_l_name.replace('ENCODE_TF_ChIP-seq_2015','...temp...').replace(
		'ChEA_2016','ENCODE_TF_ChIP-seq_2015').replace('...temp...','ChEA_2016')
	f_name = old_f_name.replace('ENCODE_TF_ChIP-seq_2015','...temp...').replace(
		'ChEA_2016','ENCODE_TF_ChIP-seq_2015').replace('...temp...','ChEA_2016')
	prefix = 'from_' + l_name + '_to_' + f_name

	#If the training data file has already been created, load it.
	fisher_classifier_fname = prefix + '_fisher_classifier_data.csv' 
	if os.path.isfile(fisher_classifier_fname): 
		print('using old file for ', fisher_classifier_fname)
		train_df = open_csv(fisher_classifier_fname)

	#Otherwise, create it. 
	else:
		print('creating ', fisher_classifier_fname)

		#Load the gmt files.
		train_l = convert_gmt('df', l_name)
		train_f = convert_gmt('df', f_name)

		#Load the files for scores which you want to train the classifier on. 
		os.chdir('..')
		os.chdir('results')
		rf_df = open_csv(prefix + '_RandomForest.csv')
		gini_df = open_csv(prefix + '_Impurity_Gini.csv')
		xgb_df = open_csv(prefix + '_XGBoost.csv')
		os.chdir('..')
		os.chdir('libs')

		#Create the training data by getting the scores for a sample of tf gene set pairs. 
		train_df = pd.DataFrame(columns=['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p','rf','gini','xgb','match'])
		#Use any score file to iterate over the tfs in the label library with at least one match in the feature library.
			#(These are the "overlaps" from get_overlaps.py .)
		for l_tf in rf_df.columns.values:
			print(l_tf)
			l_tf_genes = set(train_l.index[train_l[l_tf]])

			#Iterate over the feature library tfs. Sample half of the non-matches (for speed) and all of the matches.
			for f_tf in rf_df.index.values:
				match = clean(l_tf) == clean(f_tf)
				if rand(0,1) > .5 or match:
					f_tf_genes = set(train_f.index[train_f[f_tf]])

					#Get Fisher contingency table and p values. 
					a = len(f_tf_genes & l_tf_genes)
					b = len(f_tf_genes) - a
					c = len(l_tf_genes) - a
					d = 20000 - a - b - c
					o,p = stats.fisher_exact([[a,b],[c,d]], alternative='greater')

					#Get other scores. 
					rf = rf_df.at[f_tf,l_tf]
					gini = gini_df.at[f_tf,l_tf]
					xgb = xgb_df.at[f_tf,l_tf]

					#Append these values as a row. 
					result = pd.Series([a,b,c,d,a/d,a/(a+b),a/(a+c),p,rf,gini,xgb,match],
						index=['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p','rf','gini','xgb','match'],
						name = l_tf + ', ' + f_tf)
					train_df = train_df.append(result)

					#Since our sample is going to be very large, we will save the dataframe after every 100 pairs. 
					if train_df.shape[0] > 100:
						print('100 mark !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
						if not os.path.isfile(fisher_classifier_fname): train_df.to_csv(fisher_classifier_fname, sep='\t')
						else: 
							old = open_csv(fisher_classifier_fname)
							old = old.append(train_df)
							old.to_csv(fisher_classifier_fname, sep='\t')
						#Reset train_df .
						train_df = pd.DataFrame(columns=['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p','rf','gini','xgb','match'])

		train_df.to_csv(fisher_classifier_fname, sep='\t')

	#Use train_df to build the classifier(s).
	clf1 = XGBClassifier(random_state = 73017, n_estimators=500, learning_rate=.02)
	print('fitting classifier')
	#skip xgboost: too slow.
	clf1.fit(train_df.iloc[:,:-2], train_df['match'])
	print([i for i in zip(['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p','rf','gini','xgb','match'],clf1.feature_importances_)])
	os.chdir('..')
	os.chdir('results')

	#Return the classifier(s).
	return clf1
