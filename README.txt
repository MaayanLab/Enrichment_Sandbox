===A brief guide:===

How to use:
*Run setup.py .
*get_rankings.py has a function called get_methods(). In here, create a pandas.DataFrame with methods and parameters from enrichment_methods.py which you want to use.
*Run get_rankings.py .
*Run evaluate_rankings.py to plot the results.

Variable names:
*tf(s): transcription factor(s)

Other info:
The 'label library' is the library whose tf gene sets are being inputted for enrichment analysis.
The 'feature library' is the library which is being used as the background for enrichment analysis.
In other words, we take each tf from the label library, and perform enrichment with it using all of the feature library.
(Think of ML notation: genes are the samples, and membership within each 'feature library' gene set are the features. We want to build a classifier which can label whether or not the gene belongs to the 'label library' tf gene set. From this classifier, we extract the feature importances to get the list of enriched 'feature library' tfs.)

Results are stored as csv files for each specified method & param combination, for each library pair.
The columns of these files are the label library tfs. The indices are the feature library tfs. The cell values are the scores given by the enrichment method.
Ideally, the experiments with the best scores will correspond to the same the tf which was used to perform enrichment - this is how we benchmark the different enrichment methods!