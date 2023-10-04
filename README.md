# Sample code for "Weighted biomarker variability in joint analysis of longitudinal and time-to-event data"
Files in this repository are R codes for simulation studies in our paper "Weighted biomarker variability in joint analysis of longitudinal and time-to-event data" submitted to Annals of Applied Statistics.

Article: Weighted Biomarker Variability in Joint Analysis of Longitudinal and Time-to-Event Data

Authors: Chunyu Wang, Jiaming Shen, Christiana Charalambous, and Jianxin Pan.

This folder contains R commands and functions required to run the simulations presented in the paper. Three main R scripts are listed: 'parallel_weighted_case_1.R', 'parallel_weighted_case_3+.R', 'parallel_weighted_case_4+.R', along with their source files.

'parallel_weighted_case_1.R' along with its source file 'source_case1.R' works for the setting where both $\sigma_1$ and $\sigma_2$ are unknown and required to be estimated from the data. So 'parallel_weighted_case_1.R' and 'source_case1.R' are applicable to the estimation of case1, case2, case3 and case4 showed in the paper by changing the values of sig1ture and sig2ture in the gendat function. 

 'parallel_weighted_case_3+.R' along with its source file 'source_case3+.R' corresponds to the case 3* described in the paper: the true value of sigma2 is 10 (sig2true=10) but an approximate model without sigma2 is used to fit the dataset generated from case 3.   

 'parallel_weighted_case_4+.R' along with its source file 'source_case4+.R' corresponds to the case 4* described in the paper: the true value of sigma1 is 10 (sig1true=10) but an approximate model without sigma1 is used to fit the dataset generated from case 4. 

Due to data protection and confidentiality, the code for analysing the MRC elderly trial is not included. 
