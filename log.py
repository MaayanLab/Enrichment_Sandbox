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

RandomForest
	[m.ML_wrapper, (RandomForestClassifier, train_group, features, 73017)]

	def ML_wrapper(l_tf_genes, method, train_group, features, random_state):
	#This is a wrapper for sklearn.ensemble methods.
	target = [str(x) in l_tf_genes for x in train_group.index.values]
	clf = method(random_state = random_state)
	clf.fit(train_group[features], target)
	if method == LinearSVC: return [-abs(x) for x in clf.coef_]
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
ForestFisherCutoff.05
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
ForestFisherCutoffV2.05
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
	def InfoGain(l_tf_genes, f_matrix, info_metric):
		'''Return the tfs with descending info gain for a decision tree of length one.'''

		def Info(p,n, info_metric):
			if info_metric == 'Entropy':
				if 0 in [p,n]: return 0
				else:
					a, b = p/(p+n), n/(p+n)
					return - a * log(a,2) - b * log(b,2)

			elif info_metric == 'Gini':	
				return p/(p+n) * n/(p+n)

			else: raise Exception('Invalid info_metric.')

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

			results[column] = ((p_in + n_in) * Info(p_in, n_in, info_metric) + (p_out + n_out) * Info(p_in, n_in, info_metric)) / (p + n)
	return(list(results))
