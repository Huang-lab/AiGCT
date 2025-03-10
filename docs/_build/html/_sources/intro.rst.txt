Introduction
============

AI Genomics CollecTive (AIGCT) is a platform for systematically 
evaluating ML/AI models of variant effects across the spectrum of 
genomics-based precision medicine.

It consists of a database of variant effect data organized into categories
based on the source of the data. We use the term, task, to refer to a category.
For each task  we have the following data:

* Binary label for each variant
* VEP scores for each variant for a variety of different VEP's

The tasks available are:

* ClinVar
* Cancer - Hotspot and MSK Passenger
* ASD - Autism Spectrum Disorder
* DDD - Deciphering Developmental Disorders
* CHD - Congential Heart Disease

The platform provides a python package that supports the following functions:

* Browse and retrieve data from the database
* Evaluate the performance of a user AI/ML VEP's against
  labels and VEP's stored in our database.
* Generate a high summary report of the analysis results
* Extract the results of a performance evaluation to files for further analysis.
* Render the results in publication ready plots and figures.
