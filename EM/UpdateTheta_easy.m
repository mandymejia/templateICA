function [theta_new, subICmean, subICvar] = UpdateTheta_easy (Y, theta, C_matrix_diag, tempICmean, tempICvar)
% [theta_new, subICmean, subICvar] = UpdateTheta_easy (Y, theta, C_matrix_diag, tempICmean, tempICvar)
%
%
% %%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%
%
% Y (QxV matrix) : PCA-prewhitened fMRI data 
% theta (structure) : current parameter estimates
% theta.A           : (QxQ) mixing matrix 
% theta.nu0_sq      : (1x1) residual variance from first level
% C_matrix_diag (Qx1) : diagonal elements of matrix proportional to residual variance.  
% tempICmean (QxV matrix) : mean of each IC in template
% tempICvar (QxV matrix)  : between-subject variance of each IC in template
%
% %%%%%%%%%%%%%%% OUTPUTS %%%%%%%%%%%%%%%
%
% theta_new (structure) - updated parameter estimates
% subICmean - estimates of subject-level ICs
% subICvar - variance of subject-level ICs (for inference)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% E STEP  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First calculate the posterior probability of ICs, conditional on current parameter values

V = size(Y, 2);
Q = size(Y, 1); 

theta_new.A         = zeros(Q, Q);
theta_new.nu0_sq    = 0;

A_ProdPart1         = zeros(Q, Q);  %first part of the product format for Ai in M-step
A_ProdPart2         = zeros(Q, Q);  %second part of the product format for Ai in M-step

subICmean = zeros(Q, V); %subject-specific ICs mean
subICvar = zeros(Q, Q, V); %subject-specific ICs var

A   = theta.A ; 
C_inv = diag(1./C_matrix_diag); %QxQ
Sigma0     =   diag(C_matrix_diag.*theta.nu0_sq); %first-level covariance matrix
Sigma0_inv = diag(1./diag(Sigma0));  %first-level inverse covariance matrix


%Quantities needed for posterior probabilities of s(v)
miu_sv_part2a = (A' * C_inv * Y)./theta.nu0_sq; %mu = Cov*(part2a + part2b)
cov_sv_part1 = (A' * C_inv * A)./theta.nu0_sq;

%LOOPING OVER ALL VOXELS (COULD USE A SAMPLE TO START)
%tic
for v = 1:V

    %second-level variance of first L components
    Sigma_v = diag(tempICvar(:,v));
    Sigma_v_inv = diag(1./tempICvar(:,v));

    %E-Step: generating all moments for updating 
    Y_v  = Y(:, v);        

    %%%%% Save posterior moments of s(v)
    cov_sv = inv(cov_sv_part1 + Sigma_v_inv);
    miu_sv = cov_sv * (miu_sv_part2a(:,v) + Sigma_v_inv * tempICmean(:,v));
    miu_sv_svT = cov_sv + miu_sv * miu_sv';
    subICmean(:,v) = miu_sv;
    subICvar(:,:,v) = cov_sv;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% E STEP COMPLETE, BEGIN M STEP %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
     %%% UPDATE A
     %sum part1 and part2 over voxels, then multiply at the end to update A
     A_ProdPart1 = A_ProdPart1 + Y_v * miu_sv';
     A_ProdPart2 = A_ProdPart2 + miu_sv_svT; 
     
     %%% UPDATE nu_0^2
     theta_new.nu0_sq = theta_new.nu0_sq + (Y_v' * C_inv * Y_v) - 2*(Y_v' * C_inv * A * miu_sv) + trace(A' * C_inv * A * miu_sv_svT);

end; %FINISH LOOP OVER VOXELS
%toc

clear miu_sv miu_sv_svT miu_sv1 

%%% UPDATE PARAMETERS (M-STEP)
theta_new.nu0_sq = 1/(V*Q) * theta_new.nu0_sq;
theta_new.A = A_ProdPart1 * inv(A_ProdPart2); %new Ai
P = real(inv((theta_new.A' * theta_new.A)^(1/2))); 
theta_new.A = theta_new.A * P; % symmetric orthogonalization







