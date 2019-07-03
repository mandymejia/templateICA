function [template_mean, template_var] = estimate_templates(meas1, meas2)

  % Estimates the template mean and variance based on repeated 
  % estimates of subject-level ICs. These can be noisy as long
  % as enough subjects are included, since this function estimates 
  % and removes the noise variance.  
  % 
  % meas1 - (N x Q x V) array containing estimates of Q ICs for N subjects
  % meas2 - (N x Q x V) array containing estimates of Q ICs for N subjects
  % N = number of subjects
  % Q = number of ICs
  % V = number of locations (voxels or vertices)
  %
  % Note:
  % 1. The estimates meas1 and meas2 should be independent, i.e. from repeated 
  %    fMRI runs or from the same fMRI run, split down the middle
  % 2. The Q ICs of meas1 and meas2 must be matched.  The easiest way to 
  %    do this is to use a common set of group ICs for both sets of visits or sessions.


  %% ESTIMATE MEAN

  mean_meas1 = squeeze(nanmean(meas1,1));
  mean_meas2 = squeeze(nanmean(meas2,1));
  template_mean = (mean_meas1 + mean_meas2)/2;

  %% ESTIMATE SIGNAL (BETWEEN-SUBJECT) VARIANCE 

  % total variance
  var_tot1 = squeeze(nanvar(meas1,0,1));
  var_tot2 = squeeze(nanvar(meas2,0,1));
  var_tot = (var_tot1 + var_tot2)/2;

  % noise (within-subject) variance
  meas_diff = meas1 - meas2;
  var_noise = (1/2)*squeeze(nanvar(meas_diff,0,1));

  % signal (between-subject) variance
  template_var = var_tot - var_noise;
  template_var(template_var < 0) = 0;  







