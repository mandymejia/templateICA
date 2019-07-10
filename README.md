# templateICA
MATLAB toolbox implementing template independent component analysis

## templateICA_calc_nuisanceICs.m - helper function 
Called by Example.m

Estimates nuisance ICs a-priori using the Infomax algorithm and the ICASSO toolbox with randomized initial conditions to ensure stable ICs. Requires the [Group ICA of fMRI Toolbox](http://mialab.mrn.org/software/gift).

Uses the following files (which you should not have to modify):
-	icatb_v_pca_quiet.m - compute covariance matrix and  calculate the eigen values and eigen vectors without printing progress messages to the calling screen
-	icatb_v_whiten.m - calculate the whitening and dewhitening matrices (straight from GIFT)
-	icatb_icassoEst_quiet.m
-	icatb_icaAlgorithm_quiet.m
-	icatb_runica_quiet.m
-	getStableEStimates.m
-	changeSignOfComponents.m
