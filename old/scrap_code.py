def get_ARCHS4_correlation_matrices(lib):
	'''
	FORMERLY USED IN setup.py
	Creates .h5 files with correlation matrices between the intersection of genes from ARCHS4 and imput gmt lib files.
	libs : list
		Contains the names of the gmt files.
	'''

	new_fname = lib + '_ARCHS4_corr.h5'
	#r+ instead of w in case the operation was interrupted before, and the file was already partially created. 
	new_file = h5py.File(new_fname, 'r+')
	lib_genes = set(open_gvm(lib + '_gvm.csv').index)
	#Create a correlation matrix for both human genes and mouse genes. 
	for organism in ['human', 'mouse']:
		#In case the file was already partially created.
		if organism in list(new_file[...]): continue

		print(lib, organism)
		ARCHS4 = h5py.File(organism + '_matrix.h5', 'r')
		#Note the index contains the gene names. This is because we will use them to obtain the indices.
		ARCHS4_genes = pd.Series(range(len(ARCHS4['meta']['genes'])), index=ARCHS4['meta']['genes'])
		print(len(ARCHS4_genes))

		#Get the overlapping genes, and use them to get the overlapping indices.
		overlap_genes = {str(x).encode() for x in lib_genes} & set(ARCHS4_genes.index)
		print(len(overlap_genes))
		overlap_indices = ARCHS4_genes[overlap_genes].sort_values()

		#Get a sub-matrix by indexing only on the genes which were also in the gmt file.
		data = pd.DataFrame(ARCHS4['data']['expression'][overlap_indices,:], index=overlap_indices.index)
		print('got data')
		data = data.transpose()
		print(data.shape)

		#Remove columns with all zero values - their Pearson correlation coefficients would be undefined. 
		data = data.loc[:, (data != 0).any(axis=0)]
		print(data.shape)
		
		#Save the genes to the .h5 file.
		genes = new_file.create_dataset(organism + '/meta/genes', data = list(data.columns))
		print('got genes')

		#Obtain and save the correlation matrix to the .h5 file.
		R = data.corr()
		print('got R', R.shape)
		corr_matrix = new_file.create_dataset(organism + '/data/correlation', data = R.values)
		print('saved R')
		
		ARCHS4.close()
	new_file.close()

def download_file(url, output_fname):
	'''
	FORMERLY USED IN setup.py
	Downloads a file piecewise, and saves it to the current working directory.
	'''
	if file_exists(output_fname): return
	r = requests.get(url, stream=True)
	with open(output_fname, 'wb') as f:
		for chunk in r.iter_content(chunk_size=1024): 
			if chunk: f.write(chunk)
	return 

def pairwise_impurity_weighted(input_geneset, slib_gvm, metric):
	'''
	FORMERLY USED IN enrichment_methods.py
	Calculates impurity for each possible pair of annotations.
	An annotation's score is an adjustable function of the ranks of its top and median scores.
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
	pairs = [(i,j) for i in slib_gvm for j in slib_gvm if str(i) > str(j)]
	
	#Store the score for each pair in ```pair_results```.
	pair_results = pd.Series(index=pd.MultiIndex.from_tuples(pairs, names=['tf1','tf2']))

	#Store the score for each transcription factor. 
	#Since there are four possible final score functions, we need four ```Series```.
	individual_results1 = pd.Series(index=slib_gvm.columns)
	individual_results2 = pd.Series(index=slib_gvm.columns)
	individual_results3 = pd.Series(index=slib_gvm.columns)
	individual_results4 = pd.Series(index=slib_gvm.columns)

	#Get the set of genes corresponding to the input library tf.
	input_geneset = set(input_geneset)

	#Get the set of genes in all of the search library.
	slib_genes = set(slib_gvm.index)

	#Classify the search library genes as either in the input library tf set (p), or not (n).
	#(This is the notation from the 1975 ID3 paper.)
	p_set = slib_genes & input_geneset
	n_set = slib_genes - p_set

	#Pre-compute the sets for each tf in slib_gvm
	annot_dict = {tf:set(slib_gvm.index[slib_gvm[tf]]) for tf in slib_gvm}

	weights = [10,5,5,1]

	#For each search library tf pair, calculate the info gain resulting from the split.
	#(The split will produce four subsets whose union is the original sample, input_geneset).
	for i in set(pair_results.index.get_level_values('tf1')): #iterate over all unique tfs
		i_set = annot_dict[i]
		p_and_i, n_and_i = p_set & i_set, n_set & i_set
		p_not_i, n_not_i = p_set - p_and_i, n_set - n_and_i

		for j in pair_results[i].index:
			j_set = annot_dict[j]

			#Split the input library tf set using sets i and j.
			#Obtain the sizes of the resulting sets.
			#This is like obtaining the cell counts from a 2x2 contingency table. 
			psets = split(p_and_i, p_not_i, j_set)

			#Do the same for the genes in the input library not in the tf set. 
			nsets = split(n_and_i, n_not_i, j_set)

			#Pair the corresponding sets from p and n. 
			sets = list(zip(psets, nsets, weights))

			denom = weights[0]*(sets[0][0]+sets[0][1]) + weights[1]*(sets[1][0]+sets[1][1]) + weights[2]*(sets[2][0]+sets[2][1]) + weights[3]*(sets[3][0]+sets[3][1]) 

			#Calculate the impurity score and store its value. 
			pair_results[i,j] = sum([k[2] * (k[0] + k[1]) * I(k[0], k[1], metric) for k in sets]) / denom
			
	pair_results.sort_values(inplace=True)

	for tf in slib_gvm.columns:
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

def pairwise_impurity_nocset4(input_geneset, slib_gvm, metric):
	'''
	FORMERLY USED IN enrichment_methods.py
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
	pairs = [(i,j) for i in slib_gvm for j in slib_gvm if str(i) > str(j)]
	
	#Store the score for each pair in ```pair_results```.
	pair_results = pd.Series(index=pd.MultiIndex.from_tuples(pairs, names=['tf1','tf2']))

	#Store the score for each transcription factor. 
	#Since there are four possible final score functions, we need four ```Series```.
	individual_results1 = pd.Series(index=slib_gvm.columns)
	individual_results2 = pd.Series(index=slib_gvm.columns)
	individual_results3 = pd.Series(index=slib_gvm.columns)
	individual_results4 = pd.Series(index=slib_gvm.columns)

	#Get the set of genes corresponding to the input library tf.
	input_geneset = set(input_geneset)

	#Get the set of genes in all of the search library.
	slib_genes = set(slib_gvm.index)

	#Classify the search library genes as either in the input library tf set (p), or not (n).
	#(This is the notation from the 1975 ID3 paper.)
	p_set = slib_genes & input_geneset
	n_set = slib_genes - p_set

	#Pre-compute the sets for each tf in slib_gvm
	annot_dict = {tf:set(slib_gvm.index[slib_gvm[tf]]) for tf in slib_gvm}

	#For each search library tf pair, calculate the info gain resulting from the split.
	#(The split will produce four subsets whose union is the original sample, input_geneset).
	for i in set(pair_results.index.get_level_values('tf1')): #iterate over all unique tfs
		i_set = annot_dict[i]
		p_and_i, n_and_i = p_set & i_set, n_set & i_set
		p_not_i, n_not_i = p_set - p_and_i, n_set - n_and_i

		for j in pair_results[i].index:
			j_set = annot_dict[j]

			#Split the input library tf set using sets i and j.
			#Obtain the sizes of the resulting sets.
			#This is like obtaining the cell counts from a 2x2 contingency table. 
			psets = split(p_and_i, p_not_i, j_set)

			#Do the same for the genes in the input library not in the tf set. 
			nsets = split(n_and_i, n_not_i, j_set)

			#Pair the corresponding sets from p and n. 
			sets = zip(psets, nsets)

			#Calculate the impurity score and store its value. 
			pair_results[i,j] = sum([(k[0] + k[1]) * I(k[0], k[1], metric) for k in sets])
			
	pair_results.sort_values(inplace=True)

	for tf in slib_gvm.columns:
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