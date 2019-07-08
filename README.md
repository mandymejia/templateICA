# templateICA
MATLAB toolbox implementing template independent component analysis

## templateICA_calc_nuisanceICs.m - helper function 
Called by Example.m

Estimates nuisance ICs a-priori using the Infomax algorithm and the ICASSO toolbox with randomized initial conditions to ensure stable ICs. Requires the [Group ICA of fMRI Toolbox](http://mialab.mrn.org/software/gift).

Usage: [A_nuis, W, S_nuis, skew, iq] = templateICA_calc_nuisanceICs(ts_resid, Q2, nvox), where

Inputs:
-	ts_resid: single-subject data matrix of size TxV
-	Q2: number of nuisance ICs
-	nvox: number of voxels

Outputs:
-	A_nuis: mixing matrix for nuisance ICs
-	W: unmixing matrix
-	S_nuis: Q2 x nvox matrix of nuisance component spatial maps
-	skew: tells you whether sign of component flipped so that most voxels positively contributing to it
-	iq: numerical vector of length Q2; indicates quality of each component based on repeated runs of Infomax algorithm (usually only use components with iq>.9)

Uses the following files (which you should not have to modify):
-	icatb_v_pca_quiet.m - compute covariance matrix and  calculate the eigen values and eigen vectors without printing progress messages to the calling screen
-	icatb_v_whiten.m - calculate the whitening and dewhitening matrices (straight from GIFT)
-	icatb_icassoEst_quiet.m
-	icatb_icaAlgorithm_quiet.m
-	icatb_runica_quiet.m
-	getStableEStimates.m
-	changeSignOfComponents.m
