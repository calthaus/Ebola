###################################################################
#   R codes and data file for modeling                            #
#   sexual transmission of Ebola virus disease.                   #
#                                                                 #
#  Accompanying:                                                  #
#                                                                 #
#  Potential impact of sexual transmission on Ebola virus         #
#  epidemiology: Sierra Leone as a case study.                    #
#  JL Abbate, C-L Murall, H Richner, CL Althaus                   #
#  2016 PLoS Neglected Tropical Diseases                          #
#                                                                 #
#  PrePrint available bioRxiv http://dx.doi.org/10.1101/031880    #
#  December 11, 2015                                              #        
###################################################################

Both data files for Sierra Leone case counts are based on those downloaded from 
WHO Ebola data and statistics: Data on new cases per epi week for Sierra Leone
on 18 Nov 2015 http://apps.who.int/gho/data/node.ebola-sitrep.ebola-country-SLE-20151118


Data files and codes required to replicate our methods are found 
in this folder as such: 


1. Parameter Estimation for SEIR model 

Code: “parameter_estimation_MLE.R”
Data: “data_WHO_SierraLeone_mle.csv”


2. Fig 2

Code: “Fig2_determ_etaAnalysis.R”
Data: “data_WHO_SierraLeone_updated.csv”


3. Fig 3

Code: “Fig3_determ_alphaAnalysis.R”
Data: “data_WHO_SierraLeone_updated.csv”


4. Supporting Information S3 Fig (Deterministic Sensitivity Analysis)

Code: “determ_sensitivity.R”  (this generates results only, not final figure)


5. Fig 4 (Monte Carlo Simulations)

Code: “stoch_NoSTI.R”
Code: “stoch_STI_3mo.R”
Code: “stoch_STI_6mo.R”
Data: “data_WHO_SierraLeone_updated.csv”


6. Fig 4 (Tail Histograms)
   & Supporting Information S4 Fig (Stochastic simulations varying Eta, Violin plots)

Code: “stoch_summary_stats.R”
Data: “data_WHO_SierraLeone_updated.csv”
Data: “Stoch_noSTI_summaries_1000.csv”
Data: “Stoch_STI_3mo_summaries_1000.csv”
Data: “Stoch_STI_6mo_summaries_1000.csv”

Data: “Stoch_STI_3mo_summaries_1000_eta0005.csv”
Data: “Stoch_STI_6mo_summaries_1000_eta0005.csv”


The last two data files here were generated
using the exact same simulation files as above (for Fig 4), but changing 
eta = 0.001 to eta=0.0005 on line 49 in each of the two STI files.
