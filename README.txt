===A brief guide:===

How to use:
*Run setup.py .
*get_rankings.py has a function called get_methods(). In here, create a pandas.DataFrame with methods and parameters from enrichment_methods.py which you want to use.
*Run get_rankings.py .
*Run evaluate_rankings.py to plot the results.

Naming conventions:
*f: file
*df: pandas.DataFrame
*lib: gmt library
*l: the 'label' gmt library
*f: the 'feature' gmt library
*experiment: The column values of the original gmt files. These could be transcription factors, perterbations, etc.
*tf(s): Transcription factor(s). Note: for some gmt files, not all experiments correspond to a tf.
*funct: function - specifically, an enrichment analysis function in gsea_methods.py
*params: parameters for the function
*method: a specific function & paramaters combination which is being evaluated

What is meant by 'label' and 'feature' library?:
The 'label library' is the library whose tf gene sets are being inputted for enrichment analysis.
The 'feature library' is the library which is being used as the background for enrichment analysis.
In other words, we take each tf from the label library, and perform enrichment with it using all of the feature library.
(Think of ML notation: genes are the samples, and membership within each 'feature library' gene set are the features. We want to build a classifier which can label whether or not the gene belongs to the 'label library' tf gene set. From this classifier, we extract the feature importances to get the list of enriched 'feature library' tfs.)

How are the results stored and evaluated?:
Results are stored as csv files for each specified method & param combination, for each library pair.
The columns of these files are the label library tfs. The indices are the feature library tfs. The cell values are the scores given by the enrichment method.
Ideally, the experiments with the best scores will correspond to the same the tf which was used to perform enrichment - this is how we benchmark the different enrichment methods!
[Talk about the bridge plots]