Data Processing
- S_AdultCuriosity_Preprocessing.py takes the raw data as input and arranges the data in a csv file (matdata.csv) so that it can be fitted in the matlab model.
- S_AdultCuriosity_Processing.py must be run after S_AdultCuriosity_Preprocessing.py and loopover.m to finish up data processing before analysis. 

Bayesian Model
- loopover.m takes matdata.csv as input and fits the model looping over every sequence of every participant. The computational model mainly relies on the code for Nassar et al., 2010, Journal of Neuroscience.
- frugFunNoise.m is the Bayesian model
- seFitFrugFunNoise.m is a function that can be used to fit frugFunNoise.m to the experimental data.

Data Analysis
- StatAnalysis_S_AdultCuriosity.R takes datalp.csv and EIGdata.csv as input to run the statistical analyses reported in the paper

NOTE: all scripts are available, but the data will be shared once the paper will be accepted for pubblication