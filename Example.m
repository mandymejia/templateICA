addpath(genpath('~/matlab_toolboxes/cifti-matlab/'))
addpath(genpath('~/matlab_toolboxes/templateICA/'))
addpath(genpath('~/matlab_toolboxes/GroupICATv4.0b/')) %for GIFT

data_dir = '/path/to/timeseries/'
GICA_dir = '/path/to/groupICs/'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESTIMATE TEMPLATE MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% READ IN GROUP ICS (S0)

cd(GICA_dir)
fname = 'groupICA.dscalar.nii'; %volumetric ICs and timeseries can be used instead
groupICs = ft_read_cifti(fname);  %structure with cell 'indexmax'
V = size(groupICs.x1, 1);
Q = 25;

%put group maps in matrix form
S0 = zeros(Q, V); 
for q=1:Q
    eval(strcat('S0(q,:) = groupICs.x',num2str(q),';'));
end


%%% OBTAIN SUBJECT-LEVEL IC ESTIMATES 

fname_ts1 = 'BOLD_visit1.dtseries.nii'; %file names for voxel time courses from visit 1
fname_ts2 = 'BOLD_visit2.dtseries.nii'; %file names for voxel time courses from visit 2

maps_visit1 = zeros(S,Q,V);
maps_visit2 = zeros(S,Q,V);

%for running subjects in parallel
%parpool(12, 'IdleTimeout', Inf) %set parallel pool to never time out

N = 100; %number of subjects used to estimate templates

%parfor ii=1:N
for ii=1:N

  subject_dir = strcat('subject',num2str(ii));
  
  cd(data_dir)
  cd(subject_dir)

  dat1 = ft_read_cifti(fname_ts1); %structure with field 'dtseries'
  dat2 = ft_read_cifti(fname_ts2); %structure with field 'dtseries'

  ts1 = (dat1.dtseries)'; %TxV
  ts2 = (dat2.dtseries)'; %TxV

  %perform dual regression on each visit to get estimates
  %function below centers and scales the timeseries 

  S1 = dual_reg(ts1, S0);
  S2 = dual_reg(ts2, S0);

  maps_visit1(ii,:,:) = S1;
  maps_visit2(ii,:,:) = S2;

end


%%% ESTIMATE MEAN AND VARIANCE MAPS

[template_mean, template_var] = estimate_templates(maps_visit1, maps_visit2);
S0 = template_mean;
S0t = S0';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM TEMPLATE ICA ON A NEW SUBJECT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%read in new subject's BOLD timeseries
newsubj_dir = 'subject_new'
cd(data_dir)
cd(newsubj_dir)
dat = ft_read_cifti(fname_ts1); %structure with field 'dtseries'
ts = (dat.dtseries)'; %TxV

%perform dual regression to get first noisy estimate of template ICs
[S_DR, A_DR, ts_ctr] = dual_reg(ts, S0);

%%% ESTIMATE NUISANCE ICS WITH GIFT

%Remove estimates of template ICs
ts_resid = ts - A_DR * S_DR;

%Determine number of nuisance ICs
[~, ~, ~, ~, ~, ~, ~, Q2] = dim_reduce(ts_resid, 0);
strcat(num2str(Q2),' nuisance components') 

%Run GIFT to estimate nuisance ICs
[A_nuis, W, S_nuis, skew, iq] = icatb_calculateICA_templateICA(ts_resid, Q2, V);
sd_A = std(A_nuis); %determine scale of A
A_nuis = A_nuis * diag(1./sd_A); %rescale A
S_nuis = diag(sd_A) * S_nuis; %rescale S

%Subtract nuisance ICs from original timeseries
ts_nonuis = ts - A_nuis*S_nuis;

%Dimension-reduce data matrix using SVD 
[ts_nonuis2, H, Hinv, D_Q, U_Q, sigma2_ML, C_diag] = dim_reduce(ts_nonuis, Q); 

%Set parameter starting values
A_DR_nunuis = ts_nonuis * S0t * inv(S0 * S0t); %initialize A with DR estimate
HA = H * A_DR_nunuis; %apply dimension reduction
HA = HA * real(inv((HA' * HA)^(1/2))); %orthogonalize
sd_A = std(Hinv * HA);  %standardize scale of A (after reverse-prewhitening)
HA = HA * diag(1./sd_A);
theta0 = struct('A', HA);
theta0.nu0_sq = sigma2_ML;

%Run EM algorithm
maxiter = 100;
epsilon = .005;
[theta, S, S_var, success] = EM_easy(template_mean, template_var, ts_nonuis2, theta0, C_diag, maxiter, epsilon)
A = Hinv * theta.A; %apply reverse dimension reduction 

%S (QxV) contains the template IC estimates 
%A (TxQ) contains the template IC timeseries 









