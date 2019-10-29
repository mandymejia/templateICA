# templateICA
MATLAB toolbox implementing template Independent Component Analysis (tICA)

These files support the manuscript [Template Independent Component Analysis: Targeted and Reliable Estimation of Subject-level Brain Networks using Big Data Population Priors](https://doi.org/10.1080/01621459.2019.1679638)

If running the Fast Two-Stage EM algorithm, requires the [Group ICA of fMRI Toolbox](http://mialab.mrn.org/software/gift) to estimate nuisance ICs a-priori. 

The user is responsible for reading and writing brain image files, which can be in any format (e.g. NIFTI, CIFTI, GIFTI). The Example.m script uses CIFTI files and relies on the [cifti-matlab repository](https://github.com/Washington-University/cifti-matlab). 


## templateICA.m
Main function `templateICA`, used to fit the template ICA model using an expectation-maximization (EM) algorithm, given a set of templates.  Two EM algorithms presented in the paper are implemented: the fast two-stage EM algorithm (default) and the subspace EM algorithm.  Templates must first be estimated using the `estimate_templates` function, based on a holdout or external dataset. 

## Example.m
Example script that illustrates template estimation and template ICA model fitting.  


## EM subdirectory
Contains all functions for running the fast and subspace EM algorithms, called by the main `templateICA` function.

## GIFT_GICA subdirectory
Contains functions for estimating nuisance ICs a-priori using the Infomax algorithm and the ICASSO toolbox with randomized initial conditions to ensure stable ICs. Contains the following files, called by the main `templateICA` function:
-	icatb_v_pca.m - compute covariance matrix and calculate the eigen values and eigen vectors without printing progress messages to the calling screen
-	icatb_v_whiten.m - calculate the whitening and dewhitening matrices (straight from GIFT)
-	icatb_icassoEst.m
-	icatb_icaAlgorithm.m
-	icatb_runica.m
-	getStableEStimates.m
-	changeSignOfComponents.m
