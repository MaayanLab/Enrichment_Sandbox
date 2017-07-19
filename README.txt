===A brief guide:===

How to use:
*Run setup.py
*In get_rankings.py, create a dicts with methods and parameters you would like to use in get_methods
*Run get_rankings.py
*Run evaluate_rankings.py to plot the results

Other info:
The 'label library' is the library whose transcription factor gene sets are being inputted for enrichment analysis.
The 'feature library' is the library which is being used as the background for enrichment analysis.
In other words, we take each transcription factor from the label library, and perform enrichment on it using all of the feature library.
(Think of ML notation: genes are the samples, and membership within each feature library gene set are the features. We want to build a classifier which can detect whether or not the gene belongs to the label library transcription factor gene set. From this classifier, we extract the feature importances to get the list of enriched feature library transcription factors.)

Results are stored as csv files for each specified method & param combination, for each library pair.
Within each file, for each label library transcription factor is a ranking of the feature library transcription factors.
Ideally, the transcription factors at the top of the list will be the same as the transcription factor used to perform enrichment - this is how we benchmark the different methods!