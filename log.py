Control
	def Control(l_tf_genes, f_tfs):
		'''Return the tfs in random order.'''
		return random.sample(range(len(f_tfs)), len(f_tfs))

Fisher
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

Binomial proportion
	def BinomialProportions(l_tf_genes, f_matrix):
	'''A variation of the Fisher's exact test.'''
	p = pd.Series(index=f_matrix.columns)
	for column in f_matrix:
		f_tf_genes = set(f_matrix.index[f_matrix[column]])
		a = len(f_tf_genes & set(l_tf_genes))
		b = len(f_tf_genes)
		c = len(l_tf_genes) - a
		d = f_matrix.shape[0] - b
		#print(a,b,c,d) #diagnostics
		o,p[column] =  stats.fisher_exact([[a,b],[c,d]], alternative='greater')
	return(list(p))

RandomForest
	[m.ML_wrapper, (RandomForestClassifier, train_group, features, 73017)]

	def ML_wrapper(l_tf_genes, method, train_group, features, random_state):
		#This is a wrapper for sklearn.ensemble methods.
		target = [str(x) in l_tf_genes for x in train_group.index.values]
		clf = method(random_state = random_state)
		clf.fit(train_group[features], target)
		if method == LinearSVC: return [-abs(x) for x in clf.coef_]
		else: return [-x for x in clf.feature_importances_]

Misc ML techniques

	def ML_wrapper(l_tf_genes, method, train_group, features, random_state):
		'''This is a wrapper for sklearn.ensemble methods.'''
		target = [str(x) in l_tf_genes for x in train_group.index.values]
		clf = method(random_state = random_state)
		clf.fit(train_group[features], target)
		if method == LinearSVC: return [-abs(x) for x in clf.coef_[0]]
		else: return [-x for x in clf.feature_importances_]

CombinedFF

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
			FF = [-log(max(fi, 1e-100))*(fo + 1e-10) for (fi,fo) in zip(fisher[x], forest[x])]
			CombinedFF[x] = FF
			print(CombinedFF.shape)
		CombinedFF.to_csv(output_heading+'_CombinedFF.csv', sep='\t')
		return

	if __name__ == '__main__':
		all_libs = ['CREEDS', 'ENCODE_TF_ChIP-seq_2015', 'ChEA_2016']

		os.chdir('results')

		#Iterate over each gmt pair.
		lib_pairs = [{'l':a, 'f':b} for a in all_libs for b in all_libs if a != b]
		Parallel(n_jobs=6, verbose=0)(delayed(enrichment_wrapper)(pair)for pair in lib_pairs)

CombinedFF2

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

CombinedFF3
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
			FF = [fo if fi < .5 else fo + 1 for (fi,fo) in zip(fisher[x], forest[x])]
			CombinedFF[x] = FF
			print(CombinedFF.shape)
		CombinedFF.to_csv(output_heading+'_CombinedFF3.csv', sep='\t')
		return

	if __name__ == '__main__':
		all_libs = ['CREEDS', 'ENCODE_TF_ChIP-seq_2015', 'ChEA_2016']

		os.chdir('results')

		#Iterate over each gmt pair.
		lib_pairs = [{'l':a, 'f':b} for a in all_libs for b in all_libs if a != b]
		Parallel(n_jobs=6, verbose=0)(delayed(enrichment_wrapper)(pair)for pair in lib_pairs)

Combined,Z
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

FisherAdjusted 1-10
###
	def FisherAdjusted(l_tf_genes, f_matrix, l_lib, f_lib, ARCHS4_genes_dict):
		'''Like Fisher(), but weighs p vals by gene correlation within the intersection cell of the contingency table,
		Reward for high correlation. Also, weigh this reward by the degree of overlap with the ARCHS4 library.'''
		#Get the correlation data.
		cwd = os.getcwd()
		os.chdir('..')
		os.chdir('libs')
		ARCHS4 = h5py.File(l_lib + '_ARCHS4_corr.h5', 'r+')
		os.chdir('..')
		os.chdir(cwd)

		#Get variables which will be the same for each iteration, or only depend on the organism. 
		l_tf_genes = set(l_tf_genes)
		c_overlap_dict = {'human':{str(x).encode('utf-8') for x in l_tf_genes} & set(ARCHS4_genes_dict['human'].index),
			'mouse':{str(x).encode('utf-8') for x in l_tf_genes} & set(ARCHS4_genes_dict['mouse'].index)}
		if f_lib == 'ENCODE_TF_ChIP-seq_2015': organism_dict = {'hg19': 'human', 'mm9': 'mouse'}

		#For each tf, store the information collected in 'info' for use in the next iteration.
		info = pd.DataFrame(index=['p', 'o_frac', 'o_frac2', 'r', 'r2'], columns = f_matrix.columns)
		info.loc['o_frac',:] = 0
		info.loc['o_frac2',:] = 0
		for tf in list(f_matrix.columns):
			#Get the regular p val, just as we do in Fisher().
			f_tf_genes = set(f_matrix.index[f_matrix[tf]])
			#'a_genes' are the genes in both the feature library tf and the label library tf.
			#In other words, 'a_genes' is the intersection cell of the 2x2 contingency table. 
			a_genes = f_tf_genes & l_tf_genes
			a = len(a_genes)
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

				#Use 'ARCHS4_genes_dict' to get the appropriate list of ARCHS4 genes.
				ARCHS4_genes = ARCHS4_genes_dict[organism]

				#Get the overlap of ARCHS4 genes with the contingency table cells.
				b_overlap = {str(x).encode('utf-8') for x in f_tf_genes} & set(ARCHS4_genes.index)
				c_overlap = c_overlap_dict[organism]
				a_overlap = b_overlap & c_overlap
				a_l_o = len(a_overlap)
				b_l_o = len(b_overlap)
				c_l_o = len(c_overlap)
				l_o = len(b_overlap | c_overlap)
				
				#will edit comments below later.
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
						info.at['r',tf] = min((np.sum(r_vals) - a_l_o)/(b_l_o*c_l_o - a_l_o),.95)

						if len(a_overlap) > 1:
							info.at['o_frac2',tf] = min(a_l_o / (a),.95)
							a_overlap_indices = sorted(ARCHS4_genes[a_overlap])
							r_vals = ARCHS4[organism]['data']['correlation'][a_overlap_indices][:,a_overlap_indices]
							info.at['r2',tf] = min((np.sum(r_vals) - a_l_o)/(a_l_o*a_l_o - a_l_o),.95)

			#If the organism cannot be identified, ARCHS4 cannot be used.
			else: print('weird organism:', organism, tf)
		ARCHS4.close()

		#for tfs which did not have at least 2 genes in ARCHS4, set their 'r' val to the median. 
		#(We should not set to zero, since this would "punish" them in comparison to tfs with very low 'r' val.)
		r_grand_median = info.loc['r',:].median()
		if pd.isnull(r_grand_median): r_grand_median = 0
		info.loc['r',:].fillna(r_grand_median, inplace=True)

		r_grand_median = info.loc['r2',:].median()
		if pd.isnull(r_grand_median): r_grand_median = 0
		info.loc['r2',:].fillna(r_grand_median, inplace=True)

		#Get the adjusted p val for each tf. Lower adjusted p vals are still considered more significant.
		for tf in list(f_matrix.columns):
			info.at['p_adjusted1',tf] = math.log(info.at['p',tf]) * (1/(1-info.at['r',tf]) ** (100 * info.at['o_frac',tf]))
			info.at['p_adjusted2',tf] = math.log(info.at['p',tf]) - (25/(1-info.at['r',tf])) ** (info.at['o_frac',tf])
			info.at['p_adjusted3',tf] = math.log(info.at['p',tf]) - (1/(1-info.at['r',tf])) ** (info.at['o_frac',tf])
			info.at['p_adjusted4',tf] = math.log(info.at['p',tf]) - (100/(1-info.at['r',tf])) ** (info.at['o_frac',tf])
			info.at['p_adjusted5',tf] = (1 + math.log(info.at['p',tf])) * (1/(1-info.at['r',tf]) ** (info.at['o_frac',tf]))
			info.at['p_adjusted6',tf] = math.log(info.at['p',tf]) * (1/(1-info.at['r2',tf]) ** (100 * info.at['o_frac2',tf]))
			info.at['p_adjusted7',tf] = math.log(info.at['p',tf]) - (25/(1-info.at['r2',tf])) ** (info.at['o_frac2',tf])
			info.at['p_adjusted8',tf] = math.log(info.at['p',tf]) - (1/(1-info.at['r2',tf])) ** (info.at['o_frac2',tf])
			info.at['p_adjusted9',tf] = math.log(info.at['p',tf]) - (100/(1-info.at['r2',tf])) ** (info.at['o_frac2',tf])
			info.at['p_adjusted10',tf] = (1 + math.log(info.at['p',tf])) * (1/(1-info.at['r2',tf]) ** (info.at['o_frac2',tf]))


		return [info.loc['p_adjusted1',:], info.loc['p_adjusted2',:], info.loc['p_adjusted3',:], info.loc['p_adjusted4',:], info.loc['p_adjusted5',:], 
			info.loc['p_adjusted6',:], info.loc['p_adjusted7',:], info.loc['p_adjusted8',:], info.loc['p_adjusted9',:], info.loc['p_adjusted10',:]]

FA 11-15

	def FisherAdjusted(l_tf, l_tf_genes, f_matrix, l_lib, f_lib, ARCHS4_genes_dict):
			'''Like Fisher(), but weighs p vals by gene correlation within the intersection cell of the contingency table,
			Reward for high correlation. Also, weigh this reward by the degree of overlap with the ARCHS4 library.'''
			#Get the correlation data.
			cwd = os.getcwd()
			os.chdir('..')
			os.chdir('libs')
			ARCHS4 = h5py.File(l_lib + '_ARCHS4_corr.h5', 'r+')
			os.chdir('..')
			os.chdir(cwd)

			#Get variables which will be the same for each iteration, or only depend on the organism. 
			l_tf_genes = set(l_tf_genes)
			c_overlap_dict = {'human':{str(x).encode('utf-8') for x in l_tf_genes} & set(ARCHS4_genes_dict['human'].index),
				'mouse':{str(x).encode('utf-8') for x in l_tf_genes} & set(ARCHS4_genes_dict['mouse'].index)}
			if f_lib == 'ENCODE_TF_ChIP-seq_2015': organism_dict = {'hg19': 'human', 'mm9': 'mouse'}

			#For each tf, store the information collected in 'info' for use in the next iteration.
			info = pd.DataFrame(index=['p', 'o_frac', 'r', 'o_frac2', 'r2'], columns = f_matrix.columns)
			info.loc['o_frac',:] = 0
			info.loc['o_frac2',:] = 0
			for tf in list(f_matrix.columns):
				#Get the regular p val, just as we do in Fisher().
				f_tf_genes = set(f_matrix.index[f_matrix[tf]])
				#'a_genes' are the genes in both the feature library tf and the label library tf.
				#In other words, 'a_genes' is the intersection cell of the 2x2 contingency table. 
				a_genes = f_tf_genes & l_tf_genes
				a = len(a_genes)
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

					#Use 'ARCHS4_genes_dict' to get the appropriate list of ARCHS4 genes.
					ARCHS4_genes = ARCHS4_genes_dict[organism]

					#Get the overlap of ARCHS4 genes with the contingency table cells.
					b_overlap = {str(x).encode('utf-8') for x in f_tf_genes} & set(ARCHS4_genes.index)
					c_overlap = c_overlap_dict[organism]
					a_overlap = b_overlap & c_overlap
					a_l_o = len(a_overlap)
					b_l_o = len(b_overlap)
					c_l_o = len(c_overlap)
					l_o = len(b_overlap | c_overlap)
					
					#will edit comments below later.
					if l_o > 0: 
						#'o_frac' is the proportion of genes in the intersection cell which are also in ARCHS4.
						#Limit its value to < .95 to prevent an extreme effect when used in the adjustment formula. 
						info.at['o_frac',tf] = l_o / (b+c)
						if len(b_overlap) > 0 and len(c_overlap) > 0:
							#Get the indices of the overlapping genes, and use this to index the correlation matrix.
							b_overlap_indices = sorted(ARCHS4_genes[b_overlap])
							c_overlap_indices = sorted(ARCHS4_genes[c_overlap])
							r_vals = ARCHS4[organism]['data']['correlation'][b_overlap_indices][:,c_overlap_indices]
							#r_vals is the correlation matrix for only the overlapping tfs. 
							#Now, get the average r value for non-diagonal i.e. pairwise entries. 
							#(Each pair is duplicated across the diagonal, but this does not affect the result.)
							info.at['r',tf] = (np.sum(r_vals) - a_l_o)/(b_l_o*c_l_o - a_l_o)

							if a_l_o > 1:
								info.at['o_frac2',tf] = a_l_o / (a)
								a_overlap_indices = sorted(ARCHS4_genes[a_overlap])
								r_vals = ARCHS4[organism]['data']['correlation'][a_overlap_indices][:,a_overlap_indices]
								info.at['r2',tf] = (np.sum(r_vals) - a_l_o)/(a_l_o*a_l_o - a_l_o)

				#If the organism cannot be identified, ARCHS4 cannot be used.
				else: print('weird organism:', organism, tf)
			ARCHS4.close()

			#for tfs which did not have at least 2 genes in ARCHS4, set their 'r' val to the median. 
			#(We should not set to zero, since this would "punish" them in comparison to tfs with very low 'r' val.)

			info.to_csv('from_' + l_lib + '_to_', f_lib + '_' + l_tf.replace('/','slash') + '_fisher_adjusted_vars.csv', sep='\t')

			r_grand_median = info.loc['r',:].median()
			if pd.isnull(r_grand_median): r_grand_median = 0
			info.loc['r',:].fillna(r_grand_median, inplace=True)

			r_grand_median = info.loc['r2',:].median()
			if pd.isnull(r_grand_median): r_grand_median = 0
			info.loc['r2',:].fillna(r_grand_median, inplace=True)

			#Get the adjusted p val for each tf. Lower adjusted p vals are still considered more significant.
			for tf in list(f_matrix.columns):
				info.at['p_adjusted11',tf] = math.log(info.at['p',tf]) * (1/(1-info.at['r',tf]) ** (1 * info.at['o_frac',tf]))
				info.at['p_adjusted12',tf] = math.log(info.at['p',tf]) * ( 1 + (info.at['r',tf]) ** (1 - info.at['o_frac',tf]))
				info.at['p_adjusted13',tf] = math.log(info.at['p',tf]) * (1/(1-info.at['r2',tf]) ** (1 * info.at['o_frac2',tf]))
				info.at['p_adjusted14',tf] = math.log(info.at['p',tf]) * ( 1 + (info.at['r2',tf]) ** (1 - info.at['o_frac2',tf]))
				info.at['p_adjusted15',tf] = math.log(info.at['p',tf]) * (10/(1-info.at['r',tf]) ** (1 * info.at['o_frac',tf])) * (10/(1-info.at['r2',tf]) ** (1 * info.at['o_frac2',tf]))

			return [info.loc['p_adjusted11',:], info.loc['p_adjusted12',:], info.loc['p_adjusted13',:], info.loc['p_adjusted14',:], info.loc['p_adjusted15',:]]
			
ForestDrop
ForestDrop5
	def ML_iterative(l_tf_genes, method, it, train_group, features, random_state):
		#This is a wrapper for sklearn.ensemble methods, which chooses features recursively.
		f = list(features)
		rankings = pd.Series(index=features)
		x = 0
		while x < len(list(features)):
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

FisherForestCutoff.5
FisherForestCutoff.25
FisherForestCutoff.10
FisherForestCutoff.05
	def Fisher_ML_cutoff(l_tf_genes, method, cutoff_frac, train_group, features, random_state):
		ML_results = pd.Series(ML_wrapper(l_tf_genes, RandomForestClassifier, train_group, features, random_state), index=features)
		n_to_keep = int(len(ML_results) * cutoff_frac)
		top_features = ML_results.sort_values().index[:n_to_keep]
		p_vals_for_top_features= pd.Series(Fisher(l_tf_genes, train_group[top_features]), index=top_features)
		for f in top_features: ML_results[f] = p_vals_for_top_features[f] - 2
		return ML_results.values

ForestFisherCutoff.5
ForestFisherCutoff.25
ForestFisherCutoff.10
ForestFisherCutoff.05
	def ML_fisher_cutoff(l_tf_genes, method, cutoff_frac, train_group, features, random_state):
		fisher_results = pd.Series(Fisher(l_tf_genes, train_group), index=train_group.columns)
		n_to_keep = int(len(fisher_results) * cutoff_frac)
		new_features = fisher_results.sort_values().index[:n_to_keep]
		new_train_group = train_group[new_features]
		new_scores_for_top_features = pd.Series(ML_wrapper(l_tf_genes, RandomForestClassifier, new_train_group, new_features, random_state), 
			index=new_features)
		for f in new_features: fisher_results[f] = new_scores_for_top_features[f]
		return fisher_results.values

ForestFisherCutoffV2.5
ForestFisherCutoffV2.25
ForestFisherCutoffV2.10
ForestFisherCutoffV2.05
	def ML_fisher_cutoff_V2(l_tf_genes, method, cutoff_frac, train_group, features, random_state):
		fisher_results = pd.Series(Fisher(l_tf_genes, train_group), index=train_group.columns)
		n_to_keep = int(len(fisher_results) * cutoff_frac)
		top_features = fisher_results.sort_values().index[:n_to_keep]
		ML_scores_for_top_features = pd.Series(ML_wrapper(l_tf_genes, RandomForestClassifier, train_group, features, random_state), 
			index=features)
		for f in top_features: fisher_results[f] = ML_scores_for_top_features[f]
		return fisher_results.values


Gini
	def InfoGain(l_tf_genes, f_matrix):
		'''Return the tfs with descending Gini importances for a decision tree of length one.'''

		results = pd.Series(index=f_matrix.columns)
		l_tf_genes = set(l_tf_genes)
		f_lib_genes = set(f_matrix.index)

		p_set = f_lib_genes & l_tf_genes
		n_set = f_lib_genes - p_set
		p, n = len(p_set), len(n_set)

		for column in f_matrix:
			f_tf_genes = set(f_matrix.index[f_matrix[column]])
			p_in_set = f_tf_genes & l_tf_genes #a
			n_in_set = f_tf_genes - p_in_set #b

			p_in, n_in = len(p_in_set), len(n_in_set)
			p_out = p - p_in
			n_out = n - n_in

			results[column] = (p_in * n_in + p_out * n_out) / (p + n)
		return(list(results))

GiniImportance
	def GiniImportanceReplicate(l_tf_genes, f_matrix, info_metric):
		'''Return the tfs with descending Gini importances for a decision tree of length one.'''

		def I(p,n):
			if 0 in [p,n]: return 0
			a, b = p/(p+n), n/(p+n)
			return - a * log(a,2) - b * log(b,2)

		results = pd.Series(index=f_matrix.columns)
		l_tf_genes = set(l_tf_genes)
		f_lib_genes = set(f_matrix.index)

		p_set = f_lib_genes & l_tf_genes
		n_set = f_lib_genes - p_set
		p, n = len(p_set), len(n_set)

		for column in f_matrix:
			f_tf_genes = set(f_matrix.index[f_matrix[column]])
			p_in_set = f_tf_genes & l_tf_genes #a
			n_in_set = f_tf_genes - p_in_set #b

			p_in, n_in = len(p_in_set), len(n_in_set)
			p_out = p - p_in
			n_out = n - n_in

			results[column] =  ((p_in + n_in) * I(p_in,n_in) + (p_out + n_out) * I(p_out,n_out)) / (p + n)
		print(list(results))
		return(list(results))

InfoGain Entropy, Gini
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
			p_in corresponds to the "a" cell of the 2x2 contingency table from Fisher.
			n_in corresponds to the "b" cell of the 2x2 contingency table from Fisher.
			p_out corresponds to the "c" cell of the 2x2 contingency table from Fisher.
			n_out corresponds to the "d" cell of the 2x2 contingency table from Fisher. 
				(Although "d" is estimated as 20000 - a - b - c in our actual Fisher() function.)

			Because the starting impurity is same for each split, it is sufficient to use the weighted average
				for the two resulting nodes as the score. 
			'''
			results[column] = ((p_in + n_in) * I(p_in, n_in, metric) + (p_out + n_out) * I(p_out, n_out, metric)) / (p + n)
		return(list(results))

CombinedFG
		FF = [fo / gi for (gi,fo) in zip(gini[x], forest[x])]

ML_Fisher_features
	def get_classifier(l_name, f_name):
		os.chdir('..')
		os.chdir('libs')
		'''Under construction. Currently will not work for drug libraries.'''
		train_l_name = l_name.replace('ENCODE_TF_ChIP-seq_2015','...temp...').replace(
			'ChEA_2016','ENCODE_TF_ChIP-seq_2015').replace('...temp...','ChEA_2016')
		train_f_name = f_name.replace('ENCODE_TF_ChIP-seq_2015','...temp...').replace(
			'ChEA_2016','ENCODE_TF_ChIP-seq_2015').replace('...temp...','ChEA_2016')

		fisher_classifier_fname = 'from_' + train_l_name + '_to_' + train_f_name + '_fisher_classifier_data.csv'
		if os.path.isfile(fisher_classifier_fname): 
			print('using old file for ', fisher_classifier_fname)
			train_df = open_csv(fisher_classifier_fname)
		else:
			print('creating ', fisher_classifier_fname)
			train_l = convert_gmt('df', train_l_name)
			train_f = convert_gmt('df', train_f_name)
			train_df = pd.DataFrame(columns=['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p','match'])
			for l_tf in train_l:
				print(l_tf, train_df.shape)
				l_tf_genes = set(train_l.index[train_l[l_tf]])
				for f_tf in train_f:
					match = clean(l_tf) == clean(f_tf)
					if rand(0,1) > .95 or match:
						f_tf_genes = set(train_f.index[train_f[f_tf]])
						a = len(f_tf_genes & l_tf_genes)
						b = len(f_tf_genes) - a
						c = len(l_tf_genes) - a
						d = 20000 - a - b - c
						o,p = stats.fisher_exact([[a,b],[c,d]], alternative='greater')
						result = pd.Series([a,b,c,d,a/d,a/(a+b),a/(a+c),p,match],
							index=['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p','match'],
							name = l_tf + ', ' + f_tf)
						train_df = train_df.append(result)
			train_df.to_csv(fisher_classifier_fname, sep='\t')

		clf = RandomForestClassifier(random_state = 73017)
		print('fitting classifier')
		clf.fit(train_df.drop('match',axis=1), train_df['match'])
		print(clf.feature_importances_)
		os.chdir('..')
		os.chdir('results')
		return clf

	def ML_Fisher_features(l_tf_genes, f_matrix, classifier):
		test_df = pd.DataFrame(columns=['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p'])
		for column in f_matrix:
			f_tf_genes = set(f_matrix.index[f_matrix[column]])
			a = len(f_tf_genes & set(l_tf_genes))
			b = len(f_tf_genes) - a
			c = len(l_tf_genes) - a
			d = 20000 - a - b - c
			#print(a,b,c,d) #diagnostics
			o,p =  stats.fisher_exact([[a,b],[c,d]], alternative='greater')
			result = pd.Series([a,b,c,d,a/d,a/(a+b),a/(a+c),p],
								index=['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p'],
								name = column)
			test_df = test_df.append(result)
		return [a for (a,b) in classifier.predict_proba(test_df)]



			def get_classifiers(l_name, f_name):
				os.chdir('..')
				os.chdir('libs')
				'''Under construction. Currently will not work for drug libraries.'''
				train_l_name = l_name.replace('ENCODE_TF_ChIP-seq_2015','...temp...').replace(
					'ChEA_2016','ENCODE_TF_ChIP-seq_2015').replace('...temp...','ChEA_2016')
				train_f_name = f_name.replace('ENCODE_TF_ChIP-seq_2015','...temp...').replace(
					'ChEA_2016','ENCODE_TF_ChIP-seq_2015').replace('...temp...','ChEA_2016')

				fisher_classifier_fname = 'from_' + train_l_name + '_to_' + train_f_name + '_fisher_classifier_data.csv'
				if os.path.isfile(fisher_classifier_fname): 
					print('using old file for ', fisher_classifier_fname)
					train_df = open_csv(fisher_classifier_fname)
				else:
					print('creating ', fisher_classifier_fname)
					train_l = convert_gmt('df', train_l_name)
					train_f = convert_gmt('df', train_f_name)
					train_df = pd.DataFrame(columns=['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p','match'])
					for l_tf in train_l:
						print(l_tf, train_df.shape)
						l_tf_genes = set(train_l.index[train_l[l_tf]])
						for f_tf in train_f:
							match = clean(l_tf) == clean(f_tf)
							if rand(0,1) > .95 or match:
								f_tf_genes = set(train_f.index[train_f[f_tf]])
								a = len(f_tf_genes & l_tf_genes)
								b = len(f_tf_genes) - a
								c = len(l_tf_genes) - a
								d = 20000 - a - b - c
								o,p = stats.fisher_exact([[a,b],[c,d]], alternative='greater')
								result = pd.Series([a,b,c,d,a/d,a/(a+b),a/(a+c),p,match],
									index=['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p','match'],
									name = l_tf + ', ' + f_tf)
								train_df = train_df.append(result)
					train_df.to_csv(fisher_classifier_fname, sep='\t')

				clf1 = XGBClassifier(random_state = 73017, n_estimators=100, n_jobs=3)
				print('fitting classifier')
				clf1.fit(train_df.drop('match',axis=1), train_df['match'])
				print(clf1.feature_importances_)

				clf2 = XGBClassifier(random_state = 73017, n_estimators=500, learning_rate=.02, n_jobs=3)
				print('fitting classifier')
				clf2.fit(train_df.drop('match',axis=1), train_df['match'])
				print(clf2.feature_importances_)

				clf3 = XGBClassifier(random_state = 73017, n_estimators=100, n_jobs=3, max_depth=4)
				print('fitting classifier')
				clf3.fit(train_df.drop('match',axis=1), train_df['match'])
				print(clf3.feature_importances_)

				clf4 = XGBClassifier(random_state = 73017, n_estimators=500, learning_rate=.02, n_jobs=3, max_depth=4)
				print('fitting classifier')
				clf4.fit(train_df.drop('match',axis=1), train_df['match'])
				print(clf4.feature_importances_)

				os.chdir('..')
				os.chdir('results')
				return clf1,clf2,clf3,clf4

			df['ML_Fisher_features8'] = [m.ML_Fisher_features, (f, classifiers[0])]
			df['ML_Fisher_features9'] = [m.ML_Fisher_features, (f, classifiers[1])]
			df['ML_Fisher_features10'] = [m.ML_Fisher_features, (f, classifiers[2])]
			df['ML_Fisher_features11'] = [m.ML_Fisher_features, (f, classifiers[3])]

			def get_classifiers(l_name, f_name):
				os.chdir('..')
				os.chdir('libs')
				'''Under construction. Currently will not work for drug libraries.'''
				train_l_name = l_name.replace('ENCODE_TF_ChIP-seq_2015','...temp...').replace(
					'ChEA_2016','ENCODE_TF_ChIP-seq_2015').replace('...temp...','ChEA_2016')
				train_f_name = f_name.replace('ENCODE_TF_ChIP-seq_2015','...temp...').replace(
					'ChEA_2016','ENCODE_TF_ChIP-seq_2015').replace('...temp...','ChEA_2016')

				fisher_classifier_fname = 'from_' + train_l_name + '_to_' + train_f_name + '_fisher_classifier_data.csv'
				if os.path.isfile(fisher_classifier_fname): 
					print('using old file for ', fisher_classifier_fname)
					train_df = open_csv(fisher_classifier_fname)
				else:
					print('creating ', fisher_classifier_fname)
					train_l = convert_gmt('df', train_l_name)
					train_f = convert_gmt('df', train_f_name)
					train_df = pd.DataFrame(columns=['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p','match'])
					for l_tf in train_l:
						print(l_tf, train_df.shape)
						l_tf_genes = set(train_l.index[train_l[l_tf]])
						for f_tf in train_f:
							match = clean(l_tf) == clean(f_tf)
							if rand(0,1) > .95 or match:
								f_tf_genes = set(train_f.index[train_f[f_tf]])
								a = len(f_tf_genes & l_tf_genes)
								b = len(f_tf_genes) - a
								c = len(l_tf_genes) - a
								d = 20000 - a - b - c
								o,p = stats.fisher_exact([[a,b],[c,d]], alternative='greater')
								result = pd.Series([a,b,c,d,a/d,a/(a+b),a/(a+c),p,match],
									index=['a','b','c','d','a/d','a/(a+b)','a/(a+c)','p','match'],
									name = l_tf + ', ' + f_tf)
								train_df = train_df.append(result)
					train_df.to_csv(fisher_classifier_fname, sep='\t')

				clf1 = RandomForestClassifier(random_state = 73017, n_estimators=100)
				print('fitting classifier')
				clf1.fit(train_df.iloc[:,5:-1], train_df['match'])
				print(clf1.feature_importances_)

				clf2 = XGBClassifier(random_state = 73017)
				print('fitting classifier')
				clf2.fit(train_df.iloc[:,5:-1], train_df['match'])
				print(clf2.feature_importances_)
				os.chdir('..')
				os.chdir('results')

				df['ML_Fisher_features12'] = [m.ML_Fisher_features_2, (f, classifiers[0])]
			df['ML_Fisher_features13'] = [m.ML_Fisher_features_2, (f, classifiers[1])]

		def ML_Fisher_features_2(l_tf_genes, f_matrix, classifier):
			test_df = pd.DataFrame(columns=['a/(a+b)','a/(a+c)','p'])
			l_tf_genes = set(l_tf_genes)
			for column in f_matrix:
				f_tf_genes = set(f_matrix.index[f_matrix[column]])
				a = len(f_tf_genes & l_tf_genes)
				b = len(f_tf_genes) - a
				c = len(l_tf_genes) - a
				d = 20000 - a - b - c
				#print(a,b,c,d) #diagnostics
				o,p =  stats.fisher_exact([[a,b],[c,d]], alternative='greater')
				result = pd.Series([a/(a+b),a/(a+c),p],
									index=['a/(a+b)','a/(a+c)','p'],
									name = column)
				test_df = test_df.append(result)
			return [a for (a,b) in classifier.predict_proba(test_df)]

pairs 5,6,7,8
		if function == '1': return np.mean(scores)
		elif function == '2': return np.mean(scores[:int(len(scores)/2)])
		elif function == '3': return np.mean(scores[:int(len(scores)/10)])
		elif function == '4': return np.mean(scores[:int(len(scores)/10)]) + np.mean(scores[:int(len(scores)/2)]) * 9

		just csets 1,2,3

pairs 9,10,11,12
		if function == '1': return np.mean(scores[:int(len(scores)/10)]) #TRY THIS AGAIN?
		elif function == '2': return np.mean(scores[:int(len(scores)/50)])
		elif function == '3': return np.mean(scores[:int(len(scores)/100)])
		elif function == '4': return np.mean(scores[:int(len(scores)/1000)])

pairs 13,14,15,16
		denom = math.log((sets[0][0] + sets[0][1]) * (sets[1][0] + sets[1][1]) * (sets[2][0] + sets[0][1]) * (sets[3][0] + sets[0][1]))

		pair_results[i,j] = sum([math.log(k[0] + k[1]) * I(k[0], k[1], metric) for k in sets]) / denom

rep 1,2,3,4

			pair_results[i,j] = sum([(k[0] + k[1]) * I(k[0], k[1], metric) for k in sets])

			NO cset4

		if function == '1': return np.mean(scores[:int(len(scores)/10)]) #confirmed= replicate pair7
		elif function == '2': return np.mean(scores[:int(len(scores)/50)])
		elif function == '3': return np.mean(scores[:int(len(scores)/100)])
		elif function == '4': return np.mean(scores)

17,18,19,20

			weights = [10,5,5,1]

	#For each feature library tf pair, calculate the info gain resulting from the split.
	#(The split will produce four subsets whose union is the original sample, l_tf_genes).
	for i in set(pair_results.index.get_level_values('tf1')): #iterate over all unique tfs
		i_set = f_tf_dict[i]
		p_and_i, n_and_i = p_set & i_set, n_set & i_set
		p_not_i, n_not_i = p_set - p_and_i, n_set - n_and_i

		for j in pair_results[i].index:
			j_set = f_tf_dict[j]

			#Split the label library tf set using sets i and j.
			#Obtain the sizes of the resulting sets.
			#This is like obtaining the cell counts from a 2x2 contingency table. 
			psets = split(p_and_i, p_not_i, j_set)

			#Do the same for the genes in the label library not in the tf set. 
			nsets = split(n_and_i, n_not_i, j_set)

			#Pair the corresponding sets from p and n. 
			sets = list(zip(psets, nsets, weights))

			denom = weights[0]*(sets[0][0]+sets[0][1]) + weights[1]*(sets[1][0]+sets[1][1]) + weights[2]*(sets[2][0]+sets[2][1]) + weights[3]*(sets[3][0]+sets[3][1]) 

			#Calculate the impurity score and store its value. 
			pair_results[i,j] = sum([k[2] * (k[0] + k[1]) * I(k[0], k[1], metric) for k in sets]) / denom


more pairs
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

	#Get all possible unique pairs of transcription factors.
	pairs = [(i,j) for i in f_matrix for j in f_matrix if str(i) > str(j)]
	
	#Store the score for each pair in ```pair_results```.
	pair_results = pd.Series(index=pd.MultiIndex.from_tuples(pairs, names=['tf1','tf2']))

	#Store the score for each transcription factor. 
	#Since there are four possible final score functions, we need four ```Series```.
	individual_results1 = pd.Series(index=f_matrix.columns)
	individual_results2 = pd.Series(index=f_matrix.columns)
	individual_results3 = pd.Series(index=f_matrix.columns)
	individual_results4 = pd.Series(index=f_matrix.columns)

	#Get the set of genes corresponding to the label library tf.
	l_tf_genes = set(l_tf_genes)

	#Get the set of genes in all of the feature library.
	f_lib_genes = set(f_matrix.index)

	#Classify the feature library genes as either in the label library tf set (p), or not (n).
	#(This is the notation from the 1975 ID3 paper.)
	p_set = f_lib_genes & l_tf_genes
	n_set = f_lib_genes - p_set

	#Pre-compute the sets for each tf in f_matrix
	f_tf_dict = {tf:set(f_matrix.index[f_matrix[tf]]) for tf in f_matrix}

	weights = [10,5,5,1]

	#For each feature library tf pair, calculate the info gain resulting from the split.
	#(The split will produce four subsets whose union is the original sample, l_tf_genes).
	for i in set(pair_results.index.get_level_values('tf1')): #iterate over all unique tfs
		i_set = f_tf_dict[i]
		p_and_i, n_and_i = p_set & i_set, n_set & i_set
		p_not_i, n_not_i = p_set - p_and_i, n_set - n_and_i

		for j in pair_results[i].index:
			j_set = f_tf_dict[j]

			#Split the label library tf set using sets i and j.
			#Obtain the sizes of the resulting sets.
			#This is like obtaining the cell counts from a 2x2 contingency table. 
			psets = split(p_and_i, p_not_i, j_set)

			#Do the same for the genes in the label library not in the tf set. 
			nsets = split(n_and_i, n_not_i, j_set)

			#Pair the corresponding sets from p and n. 
			sets = list(zip(psets, nsets, weights))

			denom = weights[0]*(sets[0][0]+sets[0][1]) + weights[1]*(sets[1][0]+sets[1][1]) + weights[2]*(sets[2][0]+sets[2][1]) + weights[3]*(sets[3][0]+sets[3][1]) 

			#Calculate the impurity score and store its value. 
			pair_results[i,j] = sum([k[2] * (k[0] + k[1]) * I(k[0], k[1], metric) for k in sets]) / denom
			
	pair_results.sort_values(inplace=True)

	for tf in f_matrix.columns:
		#Get the pairs in which the tf appears first
		try:aggregated_pair_scores = list(pair_results[tf])
		except KeyError: aggregated_pair_scores = []

		#Add to that the pairs in which it appears second
		try: aggregated_pair_scores += list(pair_results[:,tf])
		except KeyError: pass

		individual_results1[tf] = final_score_function(aggregated_pair_scores, '1')
		individual_results2[tf] = final_score_function(aggregated_pair_scores, '2')
		individual_results3[tf] = final_score_function(aggregated_pair_scores, '3')
		individual_results4[tf] = final_score_function(aggregated_pair_scores, '4')

	return list(individual_results1), list(individual_results2), list(individual_results3), list(individual_results4)

def pairwise_impurity_2(l_tf_genes, f_matrix, metric):
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
		cset1 = len((c_and_i) & j_set) #i and j
		cset2 = len(c_and_i) - cset1 #i not j
		cset3 = len((c_not_i) & j_set) #j not i
		cset4 = len(c_not_i) - cset3 #neither i nor j
		return (cset1, cset2, cset3, cset4)

	def final_score_function(scores, function):
		if function == '1': return np.mean(scores[:int(len(scores)/10)])
		elif function == '2': return np.mean(scores[:int(len(scores)/50)])
		elif function == '3': return np.mean(scores[:int(len(scores)/100)])
		elif function == '4': return np.mean(scores)

	#Get all possible unique pairs of transcription factors.
	pairs = [(i,j) for i in f_matrix for j in f_matrix if str(i) > str(j)]
	
	#Store the score for each pair in ```pair_results```.
	pair_results = pd.Series(index=pd.MultiIndex.from_tuples(pairs, names=['tf1','tf2']))

	#Store the score for each transcription factor. 
	#Since there are four possible final score functions, we need four ```Series```.
	individual_results1 = pd.Series(index=f_matrix.columns)
	individual_results2 = pd.Series(index=f_matrix.columns)
	individual_results3 = pd.Series(index=f_matrix.columns)
	individual_results4 = pd.Series(index=f_matrix.columns)

	#Get the set of genes corresponding to the label library tf.
	l_tf_genes = set(l_tf_genes)

	#Get the set of genes in all of the feature library.
	f_lib_genes = set(f_matrix.index)

	#Classify the feature library genes as either in the label library tf set (p), or not (n).
	#(This is the notation from the 1975 ID3 paper.)
	p_set = f_lib_genes & l_tf_genes
	n_set = f_lib_genes - p_set

	total = f_lib_genes.shape[0]

	#Pre-compute the sets for each tf in f_matrix
	f_tf_dict = {tf:set(f_matrix.index[f_matrix[tf]]) for tf in f_matrix}

	#For each feature library tf pair, calculate the info gain resulting from the split.
	#(The split will produce four subsets whose union is the original sample, l_tf_genes).
	for i in set(pair_results.index.get_level_values('tf1')): #iterate over all unique tfs
		i_set = f_tf_dict[i]
		p_and_i, n_and_i = p_set & i_set, n_set & i_set
		p_not_i, n_not_i = p_set - p_and_i, n_set - n_and_i

		for j in pair_results[i].index:
			j_set = f_tf_dict[j]

			#Split the label library tf set using sets i and j.
			#Obtain the sizes of the resulting sets.
			#This is like obtaining the cell counts from a 2x2 contingency table. 
			psets = split(p_and_i, p_not_i, j_set)

			#Do the same for the genes in the label library not in the tf set. 
			nsets = split(n_and_i, n_not_i, j_set)

			#Pair the corresponding sets from p and n. 
			sets = zip(psets, nsets)

			#Calculate the impurity score and store its value. 

			denom = math.log((sets[0][0] + sets[0][1]) * (sets[1][0] + sets[1][1]) * (sets[2][0] + sets[0][1]) * (sets[3][0] + sets[0][1]))

			pair_results[i,j] = sum([math.log(k[0] + k[1]) * I(k[0], k[1], metric) for k in sets]) / denom
			
	pair_results.sort_values(inplace=True)

	for tf in f_matrix.columns:
		#Get the pairs in which the tf appears first
		try:aggregated_pair_scores = list(pair_results[tf])
		except KeyError: aggregated_pair_scores = []

		#Add to that the pairs in which it appears second
		try: aggregated_pair_scores += list(pair_results[:,tf])
		except KeyError: pass

		individual_results1[tf] = final_score_function(aggregated_pair_scores, '1')
		individual_results2[tf] = final_score_function(aggregated_pair_scores, '2')
		individual_results3[tf] = final_score_function(aggregated_pair_scores, '3')
		individual_results4[tf] = final_score_function(aggregated_pair_scores, '4')

	return list(individual_results1), list(individual_results2), list(individual_results3), list(individual_results4)

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
		cset1 = len((c_and_i) & j_set) #i and j
		cset2 = len(c_and_i) - cset1 #i not j
		cset3 = len((c_not_i) & j_set) #j not i
		#cset4 = len(c_not_i) - cset3 #neither i nor j
		return (cset1, cset2, cset3)

	def final_score_function(scores, function):
		if function == '1': return np.mean(scores)
		elif function == '2': return np.mean(scores[:int(len(scores)/5)])
		elif function == '3': return np.mean(scores[:int(len(scores)/10)])
		elif function == '4': return np.mean(scores[:int(len(scores)/25)])
		
	#Get all possible unique pairs of transcription factors.
	pairs = [(i,j) for i in f_matrix for j in f_matrix if str(i) > str(j)]
	
	#Store the score for each pair in ```pair_results```.
	pair_results = pd.Series(index=pd.MultiIndex.from_tuples(pairs, names=['tf1','tf2']))

	#Store the score for each transcription factor. 
	#Since there are four possible final score functions, we need four ```Series```.
	individual_results1 = pd.Series(index=f_matrix.columns)
	individual_results2 = pd.Series(index=f_matrix.columns)
	individual_results3 = pd.Series(index=f_matrix.columns)
	individual_results4 = pd.Series(index=f_matrix.columns)

	#Get the set of genes corresponding to the label library tf.
	l_tf_genes = set(l_tf_genes)

	#Get the set of genes in all of the feature library.
	f_lib_genes = set(f_matrix.index)

	#Classify the feature library genes as either in the label library tf set (p), or not (n).
	#(This is the notation from the 1975 ID3 paper.)
	p_set = f_lib_genes & l_tf_genes
	n_set = f_lib_genes - p_set

	#Pre-compute the sets for each tf in f_matrix
	f_tf_dict = {tf:set(f_matrix.index[f_matrix[tf]]) for tf in f_matrix}

	#For each feature library tf pair, calculate the info gain resulting from the split.
	#(The split will produce four subsets whose union is the original sample, l_tf_genes).
	for i in set(pair_results.index.get_level_values('tf1')): #iterate over all unique tfs
		i_set = f_tf_dict[i]
		p_and_i, n_and_i = p_set & i_set, n_set & i_set
		p_not_i, n_not_i = p_set - p_and_i, n_set - n_and_i

		for j in pair_results[i].index:
			j_set = f_tf_dict[j]

			#Split the label library tf set using sets i and j.
			#Obtain the sizes of the resulting sets.
			#This is like obtaining the cell counts from a 2x2 contingency table. 
			psets = split(p_and_i, p_not_i, j_set)

			#Do the same for the genes in the label library not in the tf set. 
			nsets = split(n_and_i, n_not_i, j_set)

			#Pair the corresponding sets from p and n. 
			sets = zip(psets, nsets)

			#Calculate the impurity score and store its value. 
			pair_results[i,j] = sum([(k[0] + k[1]) * I(k[0], k[1], metric) for k in sets])
			
	pair_results.sort_values(inplace=True)

	for tf in f_matrix.columns:
		#Get the pairs in which the tf appears first
		try:aggregated_pair_scores = list(pair_results[tf])
		except KeyError: aggregated_pair_scores = []

		#Add to that the pairs in which it appears second
		try: aggregated_pair_scores += list(pair_results[:,tf])
		except KeyError: pass

		individual_results1[tf] = final_score_function(aggregated_pair_scores, '1')
		individual_results2[tf] = final_score_function(aggregated_pair_scores, '2')
		individual_results3[tf] = final_score_function(aggregated_pair_scores, '3')
		individual_results4[tf] = final_score_function(aggregated_pair_scores, '4')

	return list(individual_results1), list(individual_results2), list(individual_results3), list(individual_results4)