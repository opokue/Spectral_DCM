
These are R codes, which create the Linear mixed model (LMM) & Function on scalar regression (FSR) analysis and bootstrapping reported in 
"Spectral Dynamic Causal Modelling of Resting-State fMRI: Relating Effective Brain Connectivity in the Default Mode Network to Genetics" 

The three files are

1. PB.R 

- This contains the r-code to run LMM Analysis for DMN4 and DMN6 PACE Data 
- and r-code to run Parametric Bootstrap for LMM. 


2. DMN4_FSR.R

- This contains the r-code for data preprocessing (data cleaning, missing value imputation, data set merging).
- PACE analysis for all edge weights on observations within 500 days.
-Function-on-scalar regression model fitting (coefficient estimation)
- F-test for FSR
-Parametric Bootstrap for FSR.



3. DMN6_FSR.R

- This contains the r-code for data preprocessing (data cleaning, missing value imputation, data set merging).
- PACE analysis for all edge weights on observations within 500 days.
-Function-on-scalar regression model fitting (coefficient estimation)
- F-test for FSR
-Parametric Bootstrap for FSR.

4. crust.R

- This contains r-code for used to obtain the clusters for the 100 SNPs.
- and also contains the code for extracting the first principal components for each cluster.







