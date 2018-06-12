import pandas as pd
import numpy as np
import os
import csv
import pickle
from scipy.sparse import coo_matrix
from itertools import compress

def open_gvm(fname):
	'''
	Returns a gvm obtained from the input file.
	fname : str
		The name of the file to be obtained as a gvm.
	'''
	if fname.endswith('csv'):
		fsize = os.path.getsize(fname)
		#For speed, read the file in chunks if it is large.
		if fsize > 1e8:
			file_chunks = pd.read_csv(fname, keep_default_na = False, sep='\t', 
			low_memory=False, encoding='Latin-1', index_col=0, chunksize=1000)
			gvm = pd.concat(file_chunks)
			gvm = gvm.replace(to_replace='',value=False).astype(bool)
		else:
			gvm = pd.read_csv(fname, keep_default_na = False, na_values=('',), sep='\t',
				low_memory=False, encoding='Latin-1', index_col=0)
			gvm = gvm.fillna(False).astype(bool)	
	elif fname.endswith('.pkl'):
			gvm = pickle.load(open(fname, 'rb'))
	else: raise ValueError('file of unknown type: ', str(fname))
	return gvm

def file_exists(fname):
	'''
	Checks if a file exists or not. Returns True or False. Prints a statement if False.
	fname : str
		The name of the file to check for existence.
	'''
	if type(fname) not in [str, bytes, os.PathLike]: return False
	if os.path.isfile(fname):
		print(fname, 'has already been created.')
		return True
	else: return False

def get_interactionlist(fname, *args):
	'''
	Returns an "interaction list" representation of a file.
	This is an nx2 pandas.DataFrame where each row gives a gene and a corresponding annotation, e.g.
		gene	annotation
		MSRA	CHEMBL406270
		TIP1	CHEMBL406270
		EPN1	CHEMBL279107
		TIP1	CHEMBL47181
		...		...
	fname : str
		The name of the file to obtain as an interaction list.
	'''

	#fivedtlibs: DrugCentral, RepurposeHub, DGIdb, DrugBank, TargetCentral
	if '1-14-18' in fname:
		interactions = pd.read_csv(fname, sep='\t', encoding='latin1')
		interactions = interactions[['TargetGeneSymbol_Entrez','gvm']]

	#STITCH
	elif '9606.protein_chemical.links.v5.0.tsv' in fname:
		CONFIDENCE_SCORES = args[0]
		print('reading file.')
		interactions = pd.read_csv('original_drug-gene_libs/9606.protein_chemical.links.v5.0.tsv', sep='\t')

		print('removing low-confidence interactions.')
		interactions = interactions.loc[interactions['combined_score'] >= min(CONFIDENCE_SCORES)]

		print('obtaining PCIDs.')
		#Ignore single vs multiple. (Previously, I confirmed that no chemical has different CIDm and CIDs combined scores.)
		interactions['chemical'] = interactions['chemical'].str.replace('CIDm', '').str.replace('CIDs', '')
		interactions['chemical'] = pd.to_numeric(interactions['chemical'])
		interactions['chemical'] = ['pc' + str(pcid) for pcid in interactions['chemical']]

		#Use lookup table to convert ENSP protein names to HGNC symbols.
		print('converting gene names.')
		gene_symbol_lookup = pd.read_csv('synonyms/ENSP_to_HGNC.csv', sep=' ')
		gene_symbol_lookup = gene_symbol_lookup.set_index('ensembl_peptide_id').to_dict()['hgnc_symbol']
		interactions['protein'] = interactions['protein'].str.replace('9606.','')
		#Keep the ENSP protein name if no HGNC symbol is available.
		interactions['protein'] = [gene_symbol_lookup.get(protein,protein) for protein in interactions['protein']]

		#Return here because we must keep the combined_score column.
		interactions = interactions.rename(columns={'chemical':'annotation', 'protein':'gene'})
		interactions = interactions.loc[~interactions['gene'].isnull().values,]
		interactions = interactions.loc[~interactions['annotation'].isnull().values,]
		interactions = interactions.drop_duplicates(subset=['gene','annotation'])
		return interactions

	##Previous version of fivedtlibs.
	#elif '_10-05-17.csv' in fname: 
	#	interactions = pd.read_csv(fname, sep=',', encoding='latin1')
	#	interactions = interactions[['TargetGeneSymbol_Entrez','DrugName']]

	##Previous version of DGIdb.
	#elif 'interactions.tsv' in fname:
	#	interactions = pd.read_csv(fname, sep='\t', encoding='latin1')
	#	interactions.loc[interactions['gene_name'].isnull().values, 'gene_name'] = interactions.loc[
	#		interactions['gene_name'].isnull().values, 'gene_claim_name']
	#	interactions = interactions[['gene_name','drug_claim_primary_name']]

	##Previous version of RepurposeHub.
	#elif 'repurposing_drugs_20170327.txt' in fname:
	#	f = pd.read_csv(fname,sep='\t',skiprows=9,encoding='latin1')
	#	f = f.loc[~f['target'].isnull().values,]
	#	interactions = []
	#	for row in f.index:
	#		geneset = f.at[row,'target'].split('|')
	#		interactions = interactions + [np.column_stack([geneset,[f.at[row,'pert_iname']] * len(geneset)])]
	#	interactions = pd.DataFrame(np.vstack(interactions))

	#DTCommons
	elif 'DTCommons_interactionlist.txt' in fname:
		interactions = pd.read_csv(fname, sep='\t')
		interactions = interactions[['gene_name', 'gvm']]

	else: raise ValueError(fname, ' was not recognized.')

	interactions.columns = ('gene','annotation')
	interactions = interactions.loc[~interactions['gene'].isnull().values,]
	interactions = interactions.loc[~interactions['annotation'].isnull().values,]
	interactions = interactions.drop_duplicates()
	return interactions

def get_genesetlist(item, item_type):
	'''
	Returns a "geneset list" representation of the input item. 
	This is a pandas.Series where the values are the annotations' geneset lists
	and the index is the annotations, e.g.
		(Index)			(Value)
		CHEMBL406270	[MSRA, TIP1, ...]
		CHEMBL279107	[ENP1, ...]
		CHEMBL47181		[TIP1, ...]
		...				...
	The input item can be a gmt file name, gvm, gvm file name, or interaction list.
	item : str (gmt file name or gvm file name), pd.DataFrame (gvm), or pd.Series (interaction list)
		The item to obtain as a geneset list.
	item_type : str
		Indicates the type of `item`. One of: "gmt_fname", "gvm", "gvm_fname", or "interactionlist"/"ilist". 
	'''
	if item_type == 'gmt_fname':
		gmt_fname = item
		with open(gmt_fname, 'r') as f:
			reader = csv.reader(f, delimiter = '\t')
			d = {row[0]:sorted(set([str(g).replace(',1.0','') for g in row[2:] if g != ''])) for row in reader}
			return pd.Series(d).sort_index()
	
	elif item_type in ['gvm_fname', 'gvm']:
		if item_type == 'gvm_fname': item, item_type = open_gvm(item), 'gvm'
		gvm = item
		genes = np.asarray(gvm.index.values)
		annotations = gvm.columns.values
		masks = [gvm[annot] == True for annot in annotations]
		d = {annot:sorted(set(genes[list(mask)])) for (annot, mask) in zip(annotations,masks)}
		return pd.Series(d).sort_index()

	elif item_type in ['interactionlist', 'ilist']: 
		return item.groupby('annotation')['gene'].apply(set).apply(sorted).sort_index()

	else: raise ValueError('unknown item type ' + item_type)

def convert_genesetlist(gslist, to, output_fname = None, verbose = False):
	'''
	Converts an input geneset list into another representation: gmt or gvm. Returns it.
	If `to == gmt` and an output file name is given, it will save the results to `output_fname`.
	If `output_fname` already exists, then the results saved to that file will be used.
	gslist : pandas.Series
		The geneset list to be converted.
	to : str
		Either 'gmt' or 'gvm'
	output_fname : str
		The name of the file to save the results to. 
	verbose : bool
		Control the frequency of print statements used when converting to gvm.

	'''
	if verbose: print('obtaining ' + output_fname)
	if to == 'gmt':
		#Create the gmt.
		gmt = [[annot] + [''] + genes for (annot,genes) in zip(gslist.index, gslist.values)]
		#Save it to the file if it does not exist yet.
		if output_fname is not None:
			if not file_exists(output_fname):
				with open(output_fname, 'w', newline='') as f:
					writer = csv.writer(f, delimiter='\t')
					for geneset in gmt: writer.writerow(geneset)
		return gmt

	elif to == 'gvm':
		#If the gvm file already exists, load it and return it.
		if file_exists(output_fname): 
			return open_gvm(output_fname)
		elif file_exists(output_fname.replace('gvm.csv','gvm.pkl')):
			return open_gvm(output_fname.replace('gvm.csv','gvm.pkl'))

		#Otherwise, create it.
		all_genes_set = {item for sublist in gslist for item in sublist}
		all_genes = pd.Series(sorted(all_genes_set))
		gslist = gslist.apply(set)
		gvm = [np.array(all_genes.isin(gs), dtype=bool) for gs in gslist]

		#Save the gvm file as a csv, or as a pickled pandas.SparseDataFrame if it is too large.
		#Transpose matrix.
		if len(gvm) < 10000:
			gvm = pd.DataFrame(gvm).transpose()
		else:
			if verbose: print('getting coo_matrix for gvm with ' + str(len(gvm)) + ' columns.')
			gvm = coo_matrix(gvm, dtype=bool).transpose()
			if verbose: print('converting coo_matrix to sparse df')
			gvm = pd.SparseDataFrame(gvm, dtype=bool, default_fill_value=False)
			if verbose: print('obtained sparse df.')
		#Format.
		gvm.index = all_genes
		gvm.columns = gslist.index
		if output_fname is not None: 
			if gvm.shape[1] < 10000: 
				gvm = gvm.replace(to_replace=False, value='')
				gvm.to_csv(output_fname, sep='\t')
			else: gvm.to_pickle(output_fname.replace('gvm.csv','gvm.pkl'))
		return gvm
	else: raise ValueError('The desired representation (`to`) is unsupported: ' + to)

def remove_small_genesets(gslist, MINIMUM_N = 5):
	'''
	Remove entries in a geneset list with less than `MINIMUM_N` elements.
	gslist : pandas.Series
		The geneset list to subset.
	MINIMUM_N : int
		The minimum number of elements for a geneset to have. 
	'''
	mask = gslist.apply(len) >= MINIMUM_N
	return gslist[mask]

def get_gmt_and_gvm(gslist, gmt_fname, gvm_fname = None):
	'''
	Shortcut which converts a geneset list to both a gmt and gvm, and saves them to files.
	gslist : pandas.Series
		The geneset list to convert.
	gmt_fname : str
		The name of the gmt file.
	gvm_fname : str
		The name of the gvm file. If unspecified, it is `gmt_fname` but with "gmt" replaced by "gvm".
	'''
	if gvm_fname is None: gvm_fname = gmt_fname.replace('gmt','gvm')
	convert_genesetlist(gslist, to='gmt', output_fname=gmt_fname)
	convert_genesetlist(gslist, to='gvm', output_fname=gvm_fname)

def combine_genesetlists(gslist1, gslist2, merge_type = 'union'):
	'''
	Merges two geneset lists with the same annotation. Returns the new geneset list.
	gslist1 : pandas.Series
	gslist2 : pandas.Series
	merge_type : str
		"union" or "intersection"
	'''
	if not len(gslist1) == len(gslist2): raise ValueError('gslist1 and gslist2 must be the same length.')
	if not all(gslist1.index == gslist2.index): 
		gslist1, gslist2 = gslist1.sort_index(), gslist2.sort_index()
		if not all(gslist.index == gslist2.index): raise ValueError(
			'gslist1 and gslist2 must have the same annotations.')

	if merge_type == 'union':
		return pd.Series(
			{i:sorted(set(e1) | set(e2)) for (i, e1, e2) in zip(gslist1.index, gslist1, gslist2)}
		).sort_index()

	elif merge_type == 'intersection':
		return pd.Series(
			{i:sorted(set(e1) & set(e2)) for (i, e1, e2) in zip(gslist1.index, gslist1, gslist2)}
		).sort_index()

	else: raise ValueError('merge_type must be union or intersection.')

def expand_gvm(gvm, expansion, output_fname):
	'''
	Expands the genesets in a gvm using the expansion library. 
	An annotation's new geneset will be the union of each of its genes' "expansion set," 
	the gene's most co-expressed or interacting genes.
	(If a gene is not in the expansion library, its "expansion set" is just itself.)
	gvm : pandas.DataFrame
		The gvm to expand.
	expansion : pandas.DataFrame
		The expansion library. Each row is indexed by the gene and 
		the column values are the most co-expressed/interacting genes.
	output_fname : str
		The name of the file to save the expanded gvm to.
	'''

	#Exit if the gvm has already been expanded.
	if(os.path.isfile(output_fname)):
		print(output_fname, 'already created.')
		return(None)

	#We will only need expansion sets for genes which are in the gvm.
	expansion = expansion[[g for g in expansion.columns if g in gvm.index]]
	expansion_genes = set(expansion.columns)
	gene_union = sorted(set(gvm.index.values).union({i for s in expansion.values for i in s if i == i}))
	
	#Pre-allocate the expansion gvm and fill annotation by annotation.
	expanded_gvm = pd.DataFrame(False, index=gene_union, columns=gvm.columns)
	for annot in expanded_gvm.columns:
		#Get this annotation's geneset.
		genes = set(compress(gvm.index.values, gvm[annot].values))
		#For each gene in this geneset, get its expansion set (most co-expressed or co-interacting).
		##Take the union of these expansion sets.
		genes_in_exp = list(genes.intersection(expansion_genes))
		if len(genes_in_exp) == 0:
			expanded_gvm.loc[genes, annot] = True
			continue
		exp_genes = {i for s in expansion[genes_in_exp].values for i in s if i == i}
		#Add the annotation's original geneset to this list.
		genes_to_add = genes.union(exp_genes)
		#Save to the gvm.
		expanded_gvm.loc[genes_to_add, annot] = True

	#Save the gvm file as a csv, or as a pickled pandas.SparseDataFrame if it is too large.
	if expanded_gvm.shape[1] < 10000: 
		expanded_gvm = expanded_gvm.replace(to_replace=False, value='')
		expanded_gvm.to_csv(output_fname, sep='\t')
	else: pd.SparseDataFrame(expanded_gvm, default_fill_value=False).to_pickle(output_fname.replace('gvm.csv','gvm.pkl'))

def relabel_duplicates(l):
	'''
	Replicates pandas re-labeling of duplicate columns. For example:
	['ABC','DEF','ABC','GHI','GHI','ABC'] --> 
	['ABC.1','DEF','ABC.2','GHI.1','GHI.2','ABC.3']
	l : list of str
		The list to re-label.
	'''
	def relabel(old_label, count):
		return old_label + '.' + str(count)

	l = list(l)
	counts = {}
	for elem in l:
		if elem in counts: counts[elem] += 1
		else: counts[elem] = 1
	duplicates = [k for k in counts.keys() if counts[k] > 1]

	counts = {}
	for i in range(len(l)):
		elem = l[i]
		if elem in duplicates:
			if elem in counts: counts[elem] += 1
			else: counts[elem] = 1
			l[i] = relabel(elem, counts[elem])
			
	return(l)