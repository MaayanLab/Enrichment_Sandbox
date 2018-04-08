import pandas as pd
import numpy as np
import os
import csv

def open_gvm(fname):
	#Open the file.
	gvm = pd.read_csv(fname, keep_default_na = False, na_values=('',), sep='\t', 
		low_memory=False, encoding='Latin-1', index_col=0)
	
	##Workaroud from bug which raised error if keep_default_na=False and index_col=0 are both used.
	#gvm = pd.read_csv(fname, keep_default_na = False, na_values=('',), sep='\t', low_memory=False, encoding='Latin-1')
	#lib_name = fname.partition('_gvm.csv')[0]
	#gvm.set_index(gvm[lib_name], inplace=True)
	#gvm.drop([lib_name], axis=1, inplace=True)

	#Convert blank cells to False, and ensure bool type. 
	gvm = gvm.fillna(False).astype(bool)
	return gvm

def file_exists(fname):
	if type(fname) not in [str, bytes, os.PathLike]: return False
	'''Checks if a file exists in the directory, printing a statement if so.'''
	if os.path.isfile(fname):
		print(fname, 'has already been created.')
		return True
	else: return False

def get_interactionlist(fname, *args):
	#fivedtlibs
	if '1-14-18' in fname:
		interactions = pd.read_csv(fname, sep=',', encoding='latin1')
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

	elif '_10-05-17.csv' in fname: 
		interactions = pd.read_csv(fname, sep=',', encoding='latin1')
		interactions = interactions[['TargetGeneSymbol_Entrez','DrugName']]

	elif 'interactions.tsv' in fname:
		interactions = pd.read_csv(fname, sep='\t', encoding='latin1')
		interactions.loc[interactions['gene_name'].isnull().values, 'gene_name'] = interactions.loc[
			interactions['gene_name'].isnull().values, 'gene_claim_name']
		interactions = interactions[['gene_name','drug_claim_primary_name']]

	elif 'repurposing_drugs_20170327.txt' in fname:
		f = pd.read_csv(fname,sep='\t',skiprows=9,encoding='latin1')
		f = f.loc[~f['target'].isnull().values,]
		interactions = []
		for row in f.index:
			geneset = f.at[row,'target'].split('|')
			interactions = interactions + [np.column_stack([geneset,[f.at[row,'pert_iname']] * len(geneset)])]
		interactions = pd.DataFrame(np.vstack(interactions))

	elif 'dtc_interactionlist.txt' in fname:
		interactions = pd.read_csv(fname, sep='\t')
		interactions = interactions[['gene_name', 'gvm']]

	else: raise ValueError(fname, ' was not recognized.')

	interactions.columns = ('gene','annotation')
	interactions = interactions.loc[~interactions['gene'].isnull().values,]
	interactions = interactions.loc[~interactions['annotation'].isnull().values,]
	interactions = interactions.drop_duplicates()
	return interactions

def get_genesetlist(item, item_type):
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
		d = {annot:sorted(set(genes[mask])) for (annot, mask) in zip(annotations,masks)}
		return pd.Series(d).sort_index()

	elif item_type in ['interactionlist', 'ilist']: 
		return item.groupby('annotation')['gene'].apply(set).apply(sorted).sort_index()

	else: raise ValueError('unknown item type ' + item_type)

def convert_genesetlist(gslist, to, output_fname = None, verbose = False):
	if to == 'gmt':
		#Create the gmt.
		gmt = [[annot] + [''] + genes for (annot,genes) in zip(gslist.index, gslist.values) if len(genes) > 4]
		#Save it to the file if it does not exist yet.
		if output_fname is not None:
			if not file_exists(output_fname):
				with open(output_fname, 'w', newline='') as csvfile:
					writer = csv.writer(csvfile, delimiter='\t')
					for geneset in gmt: writer.writerow(geneset)
		return gmt

	if to == 'gvm':
		#If the gvm file already exists, load it and return it.
		if file_exists(output_fname): return open_gvm(output_fname)
		#Otherwise, create it, save it to the file, and return it.
		all_genes = sorted({item for sublist in gslist for item in sublist})
		gvm = pd.DataFrame([[gene in s for gene in all_genes] for s in gslist]).transpose()
		gvm.index = all_genes
		gvm.columns = gslist.index #These are the annotations.
		gvm = gvm.replace(to_replace=False, value='')
		if output_fname is not None: gvm.to_csv(output_fname, sep='\t')
		return gvm

def remove_few_targets(gslist, MINIMUM_N_TARGETS = 5):
	mask = gslist.apply(len) >= MINIMUM_N_TARGETS
	return gslist[mask]

def combine_genesetlists(gslist1, gslist2, merge_type = 'union'):
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