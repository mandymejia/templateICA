function [template_mean, template_var] = estimate_templates(meas1, meas2)

  % Estimates the template mean and variance based on repeated 
  % estimates of subject-level ICs. These can be noisy as long
  % as enough subjects are included, since this function estimates 
  % and removes the noise variance.
  % 
  % meas1 - (N x V) array containing IC estimates for N subjects
  % meas2 - (N x V) array containing IC estimates for N subjects
  % meas1 and meas2 should be independent, i.e. from repeated fMRI runs
  % or from the same fMRI run, split down the middle
  % N = number of subjects
  % V = number of locations (voxels or vertices)


  %% ESTIMATE MEAN

  mean_meas1 = squeeze(nanmean(meas1));
  mean_meas2 = squeeze(nanmean(meas2));
  template_mean = (mean_meas1 + mean_meas2)/2;

  %% ESTIMATE SIGNAL (BETWEEN-SUBJECT) VARIANCE 

  % Estimate total variance
  var_tot1 = squeeze(nanvar(meas1));
  var_tot2 = squeeze(nanvar(meas2));
  var_tot = (var_tot1 + var_tot2)/2;

  % Estimate noise (within-subject) variance
  meas_diff = meas1 - meas2;
  var_noise = (1/2)*squeeze(nanvar(meas_diff));

  % Estimate signal variance
  template_var = var_tot - var_noise;
  template_var(template_var < 0) = 0;  







