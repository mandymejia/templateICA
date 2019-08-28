function [S, A, Q_nuis] = templateICA(dat, tempICmean, tempICvar, flag, maxQ, maxiter, epsilon)
% Usage: [S, A, Q_nuis] = templateICA(dat, tempICmean, tempICvar, flag, maxQ, maxiter, epsilon)
%
% Performs centering, dimension reduction, gets starting values, and calls EM_easy or EM_subspace
%
% 
% ARGUMENTS
%
% dat        - (TxV) matrix of fMRI data
% tempICmean - (LxV) matrix of template means
% tempICvar  - (LxV) matrix of template variances
% flag       - If flag=1, use subspace EM, otherwise use fast EM
% maxQ       - Maximum number of ICs (template+nuisance) to identify (L <= maxQ <= T)
% maxiter    - maximum number of EM iterations
% epsilon    - smallest proportion change between iterations (e.g. .01 or 1%)
%
%
% RETURNS
%
% S      - (QxV) matrix of IC estimates (template|nuisance)
% A      - (TxQ) mixing matrix
% Q_nuis - how many nuisance ICs estimated

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHECK AND SET ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if extra arguments given
if(nargin > 7) 
  error('Too many arguments supplied');
end;

%check that first three arguments are legit
stuff_missing = ~exist('dat','var') || ~exist('tempICmean','var') || ~exist('tempICvar','var');
if stuff_missing
  error('Must supply data, template mean and template variance arguments');
else 
  stuff_empty = isempty(dat) || isempty(tempICmean) || isempty(tempICvar);
  stuff_notmatrix = ~ismatrix(dat) || ~ismatrix(tempICmean) || ~ismatrix(tempICvar);
  if stuff_empty || stuff_notmatrix
    error('Please give me non-empty matrices for data, template mean and template variance arguments');
  end
end

L = size(tempICmean,1); %number of template ICs
V = size(tempICmean,2); %number of voxels/vertices
T = size(dat,1);

if(T > V) error('More time points than voxels, check matrix orientation'); end
if(L > V) error('More ICs than voxels, check matrix orientation'); end
if(L > T) error('More ICs than time points, check matrix orientation'); end
if(V ~= size(dat, 2)) error('The number of voxels in data and template must match'); end
if(L ~= size(tempICvar,1)) error('The number of ICs in the template mean and variance must match'); end
if(V ~= size(tempICvar,2)) error('The number of voxels in the template mean and variance must match'); end

%set flag=0 to use fast EM algorithm unless flag=1
if ~exist('flag','var') || isempty(flag)
  flag=0;
end

if(flag == 1)
  disp('Using subspace EM algorithm')
else 
  disp('Using fast EM algorithm')
end

%set maxQ to T if not specified or is more than T
if ~exist('maxQ','var') || isempty(maxQ)
  maxQ = T;
elseif maxQ > T
  maxQ = T;
end

%check that maxQ makes sense
if(maxQ < L) error('maxQ must be at least the number of template ICs'); end
if(maxQ == L) warning('Since maxQ is the number of template ICs, no nuisance ICs will be estimated. Template IC estimates may be contaminated as a result.'); end;

%set maxiter and epsilon if not supplied
if (~exist('maxiter','var') || isempty(maxiter)) 
  maxiter = 100; 
  disp('You did not provide a value for maxiter, so I am setting it to 100.')
end;
if (~exist('epsilon','var') || isempty(epsilon)) 
  epsilon = 0.01; 
  disp('You did not provide a value for epsilon, so I am setting it to 0.01.')
end;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ESTIMATE AND DEAL WITH NUISANCE ICS (unless maxQ = L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Versions of the fMRI data below:
%dat – original data, becomes centered data through dual_reg function
%dat2 - data after removal of template ICs, used to estimate nuisance ICs
%dat3 - data used in EM algorithm
%     - for fast EM: data after removal of nuisance ICs
%     - for subspace EM: same as dat
%dat4 - dat3 after dimension reduction, used in EM algorithm

if(maxQ > L)

  %Get initial estimate of A and S associated with template ICs using dual regression
  [S_temp_DR, A_temp_DR, dat] = dual_reg(dat, tempICmean);

  %To estimate number of nuisance ICs, first estimate and remove template ICs (more accurate than estimating the total number of ICs)

  fprintf('\nESTIMATING NUMBER OF OF NUISANCE COMPONENTS...:')

  %Remove A*S from data to estimate nuisance ICs
  dat2 = dat - A_temp_DR * S_temp_DR; %dat2 should contain only nuisance ICs
  [~, ~, ~, ~, ~, ~, ~, Q_nuis] = dim_reduce(dat2, 0); %estimate number of nuisance ICs
  if(Q_nuis+L > maxQ) Q_nuis = maxQ - L; end

  fprintf('%4.0f \n', Q_nuis)
  
  fprintf('\nRUNNING GIFT TO OBTAIN INITIAL ESTIMATES OF NUISANCE COMPONENTS...\n')

  %Run Infomax on dat2 to estimate nuisance ICs
  tic
  [A_nuis, ~, S_nuis, ~, ~] = icatb_calc_nuisanceICs(dat2, Q_nuis, V);
  toc
  sd_A = std(A_nuis); %determine scale of A
  A_nuis = A_nuis * diag(1./sd_A); %rescale A
  S_nuis = diag(sd_A) * S_nuis; %rescale S


  %%% DEAL WITH NUISANCE ICS (*how* depends on algorithm used)

  %%% Fast EM: Subtract nuisance components 

  if(flag~=1) 
    
    dat3 = dat - A_nuis*S_nuis; %run EM on data with nuisance ICs removed
    Q_EM = L; %number of ICs to estimate with EM (= number of PCs to keep)

  %%% Subspace EM: Use estimated nuisance ICs to get MoG starting values

  else 
    dat3 = dat; %run EM on full data

    M = 2; %number of MoG components to use (M=2 only computationally feasible option in practice)
    [miu0, sigma_sq0, pi0] = deal(zeros(M, Q_nuis)); % after reshaping: [miu_11, ..., miu_1M, miu_21, ..., miu_2M, ..., miu_L1, ..., miu_LM]'
    for(k = 1:Q_nuis)
      mog_k = MoGfit(S_nuis(k,:),M);
      %in subspace EM, last MoG component is assumed to represent signal
      if(mog_k.pi(1) < mog_k.pi(2)) %if true, first MoG component is signal, need to switch since we assume the first MoG component is background (majority of voxels)
        o = [2,1]; %switch order of MoG components
      else 
        o = 1:2; %keep order of MoG components
      end
      miu0(:,k) = mog_k.miu(o);
      sigma_sq0(:,k) = mog_k.sigma_sq(o);
      pi0(:,k) = mog_k.pi(o);    
    end
    Q_EM = L + Q_nuis; %number of ICs to estimate with EM (= number of PCs to keep)
  end

%If maxQ=L, no need to deal with nuisance ICs
else
  dat3 = dat; %we didn't need to remove any nuisance ICs, so just run EM on the original data
  flag = 0; %use fast EM since it is the same as subspace EM when there are no nuisance ICs
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DIMENSION REDUCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nPERFORMING DIMENSION REDUCTION\n')

%Y_ij2 is dimension-reduced data (QxV)
%H is dimension reduction matrix
%Hinv is reverse dimension reduction matrix
%sigma2_ML is the initial residual variance estimate, based on the averaged non-retained eigenvalues
%C_diag is the residual covariance structure induced by dimension reduction
[dat4, H, Hinv, ~, ~, sigma2_ML, C_diag] = dim_reduce(dat3, Q_EM); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INITIALIZE PARAMETER STARTING VALUES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nSETTING INITIAL PARAMETER VALUES\n')

% MIXING MATRIX 
[~, A_temp_init, ~] = dual_reg(dat3, tempICmean); %initialize A with DR estimate
HA = H * A_temp_init; %apply dimension reduction
HA = HA * real(inv((HA' * HA)^(1/2))); %orthogonalize
sd_A = std(Hinv * HA);  %standardize scale of A (after reverse-prewhitening)
HA = HA * diag(1./sd_A);
theta0 = struct('A', HA);

% RESIDUAL VAR
theta0.nu0_sq = sigma2_ML;

% MOG PARAMETERS (FOR SUBSPACE EM ALGORITHM ONLY)
if(flag == 1)
  HA_nuis = H * A_nuis;
  HA_nuis = HA_nuis * real(inv((HA_nuis' * HA_nuis)^(1/2))); %orthogonalize
  sd_A = std(Hinv * HA_nuis);  %standardize scale of A (after reverse-prewhitening)
  HA_nuis = HA_nuis * diag(1./sd_A);

  %theta0.A = [theta0.A, normrnd(0,1,[Q_EM,Q_nuis])]; %use A_nuis here?
  theta0.A = [theta0.A, HA_nuis];
  theta0.miu = reshape(miu0, [M*Q_nuis, 1]);
  theta0.sigma_sq = reshape(sigma_sq0, [M*Q_nuis, 1]);
  theta0.pi = reshape(pi0, [M*Q_nuis, 1]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% RUN EM ALGORITHM UNTIL CONVERGENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nRUNNING EM-ALGORITHM...\n')

% FAST EM ALGORITHM

if(flag ~= 1)

  %theta is the EM parameter estimates
  %S_temp is the estimated template ICs
  %S_temp_var is the estimation variance (squared standard errors) of the ICs
  %success is a flag for algorithm convergence
  tic
  [theta, S_temp, S_temp_var, success] = EM_easy(tempICmean, tempICvar, dat4, theta0, C_diag, maxiter, epsilon);
  toc
  A_temp = Hinv * theta.A; %reverse dimension reduction to get TxQ mixing matrix

  fprintf('\nRE-ESTIMATING NUISANCE COMPONENTS WITH GIFT...\n')

  %Re-estimate S_nuis and A_nuis and scale
  dat2 = dat - A_temp * S_temp;
  tic
  [A_nuis, ~, S_nuis, ~, ~] = icatb_calc_nuisanceICs(dat2, Q_nuis, V);
  toc
  sd_A = std(A_nuis); %determine scale of A
  A_nuis = A_nuis * diag(1./sd_A); %rescale A
  S_nuis = diag(sd_A) * S_nuis; %rescale S

  %Combine template and nuisance ICs
  S = [S_temp; S_nuis];
  A = [A_temp, A_nuis];

else

% SUBSPACE EM ALGORITHM

  %theta is the EM parameter estimates
  %S is the estimated template and nuisance ICs
  %S_temp_var is the estimation variance (squared standard errors) of the ICs
  %success is a flag for algorithm convergence
  tic
  [theta, S, S_var, success] = EM_subspace(tempICmean, tempICvar, dat4, theta0, C_diag, maxiter, epsilon);
  toc
  A = Hinv * theta.A; %reverse dimension reduction to get TxQ mixing matrix

end

fprintf('\nDONE!\n')

