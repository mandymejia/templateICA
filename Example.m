addpath(genpath('~/matlab_toolboxes/cifti-matlab/'))
addpath(genpath('~/matlab_toolboxes/templateICA/'))

data_dir = '/path/to/timeseries/'
GICA_dir = '/path/to/groupICs/'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESTIMATE TEMPLATE (MEAN AND VARIANCE) MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OBTAIN SUBJECT-LEVEL IC ESTIMATES 
S = 100; %number of subjects to use in template estimation

fname_ts1 = 'BOLD_visit1.dtseries.nii'; %file names for voxel time courses from visit 1; saved in subject-specific directories
fname_ts2 = 'BOLD_visit2.dtseries.nii'; %file names for voxel time courses from visit 2; saved in subject-specific directories

%initialize subject-level maps
maps_visit1 = zeros(S,Q,V);
maps_visit2 = zeros(S,Q,V);

%for running subjects in parallel
%parpool(12, 'IdleTimeout', Inf) %set parallel pool to never time out

%parfor ii=1:S
for ii=1:S
  fprintf('working on %d of %d', ii, S)

  subject_dir = strcat('subject',num2str(ii));

  dat1 = ft_read_cifti(fullfile(data_dir, subject_dir, fname_ts1)); %structure with field 'dtseries'
  dat2 = ft_read_cifti(fullfile(data_dir, subject_dir, fname_ts2)); %structure with field 'dtseries'

  ts1 = (dat1.dtseries)'; %TxV
  ts2 = (dat2.dtseries)'; %TxV

  %perform dual regression on each visit to get estimates
  %function below centers and scales the timeseries 

  S1 = dual_reg(ts1, S0);
  S2 = dual_reg(ts2, S0);

  maps_visit1(ii,:,:) = S1;
  maps_visit2(ii,:,:) = S2;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULATE MEAN AND SIGNAL VARIANCE 

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

%perform dual regression (for comparison) and center data
[S_DR, A_DR, ts_ctr] = dual_reg(ts, S0);

%perform template ICA using fast EM algorithm
[S_tICA, A_tICA, Q2] = templateICA(ts_ctr, template_mean, template_var, 0, 200); %max 200 ICs

%perform template ICA using subspace EM algorithm (WARNING: may be slow!)
[S_tICA_subspace, A_tICA_subspace, Q2] = templateICA(ts_ctr, template_mean, template_var, 1, 200); %max 200 ICs


%S (QxV) contains the template IC estimates 
%A (TxQ) contains the template IC timeseries 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WRITE OUT ESTIMATES AS CIFTI TIMESERIES FILE (USE CONNECTOME WORKBENCH TO VISUALIZE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%https://github.com/Washington-University/HCPpipelines/tree/master/global/matlab
addpath(genpath('~/matlab_toolboxes/gifti-1.6/'))
addpath('~/matlab_toolboxes/cifti/') %ciftiopen, ciftisave, ciftisavereset
wb_cmd = '~/workbench/bin_rh_linux64/wb_command';

cd(data_dir)
cd(newsubj_dir)
S_cifti = ciftiopen(fname_ts1);
S_cifti.cdata = S_DR';
ciftisavereset(S_cifti, 'S_DR.dscalar.nii', wb_cmd)

S_cifti.cdata = S_tICA';
ciftisavereset(S_cifti, 'S_tICA.dscalar.nii', wb_cmd)

S_cifti.cdata = S_tICA_subspace';
ciftisavereset(S_cifti, 'S_tICA_subspace.dscalar.nii', wb_cmd)



