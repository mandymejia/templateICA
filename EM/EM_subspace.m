function [theta, subICmean, subICvar, success, zero_post_z] = EM_subspace (tempICmean, tempICvar, Y, theta0, C_matrix_diag, maxiter, epsilon)
% [theta, subICmean, subICvar, success, zero_post_z] = TemplateICA_EM (tempICmean, tempICvar, Y, theta0, C_matrix_diag, maxiter, epsilon)
%
% EM algorith for single-subject template ICA model (easy version)
%
% Q - number of ICs in template
% V - number of voxels/vertices in data
% T - number of original fMRI time points
%
% %%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%
%
% tempICmean (QxV matrix) - mean of each IC in template
% tempICvar (QxV matrix) - between-subject variance of each IC in template
%
% Y (QxV matrix) - PCA-prewhitened fMRI data 
%
% theta0 (structure) - initial guess at parameter values
% theta.A             : (QxQ) mixing matrix 
% theta.nu0_sq        : (1x1) residual variance from first level
% theta.miu         : (M(Q-L)x1) miu_z (MoG means)
% theta.sigma_sq    : (M(Q-L)x1) sigma^2 (MoG variances)
% theta.pi          : (M(Q-L)x1) pi (MoG component probabilities)
%
% C_matrix_diag (Qx1) - diagonal elements of matrix proportional to residual variance.  
% If original fMRI timeseries has covariance sig^2*I, the prewhitened timeseries 
% achieved by premultiplying by (QxT) matrix H from PCA has diagonal covariance sig^2*HH', 
% so C_matrix_diag = diag(HH').
%
% maxiter - maximum number of EM iterations
% epsilon - smallest percent change between iterations (e.g. .001)
%
% %%%%%%%%%%%%%%% OUTPUTS %%%%%%%%%%%%%%%
%
% theta (structure) - final parameter estimates
% subICmean - estimates of subject-level ICs
% subICvar - variance of subject-level ICs (for inference)
% success - flag indicating convergence (1) or not (0)
% zero_post_z - number of voxels for which posterior probability of z is zero within subspace

Q = size(Y, 1);
V = size(Y, 2);
itr = 1;
err= 1000; %set large initial value for difference between iterations
theta = theta0;
success = 1;
tempICvar(tempICvar < .00001) = .00001; %to prevent problems when inverting covariance
zero_post_z = [];

while err > epsilon
    [theta_new, z_mode, subICmean, subICvar, error, zero_post_z_iter] = UpdateTheta_subspace(Y, theta, C_matrix_diag, tempICmean, tempICvar);
    %subICmean, subICvar are mean and variance of subject-specific ICs
    %z-mode indicates which MoG component each voxel belongs to with highest prob

    zero_post_z = [zero_post_z, zero_post_z_iter];
    if(error == 1)
        success = 0;
        disp('Fail to converge due to calculate of p(z|y)!');
        return;
    end;

    A_vec = reshape(theta.A, [Q*Q,1]);
    A_vec_new = reshape(theta_new.A, [Q*Q,1]);
    err_A = norm(A_vec_new - A_vec)/norm(A_vec);
    err_nu = abs((theta_new.nu0_sq - theta.nu0_sq))/theta.nu0_sq;
    err_miu = norm(theta_new.miu - theta.miu)/norm(theta.miu);
    err_sigma_sq = norm(theta_new.sigma_sq - theta.sigma_sq)/norm(theta.sigma_sq);
    err_pi = norm(theta_new.pi - theta.pi)/norm(theta.pi);
    err_all = [err_A, err_nu, err_miu, err_sigma_sq, err_pi]
    err = max(err_all);

    fprintf('iteration %6.0f and the difference is  %6.6f for theta \n', itr, err);  
    clear A_vec A_vec_new err_A err_nu;
    theta = theta_new;
    itr = itr + 1;
    if ( itr > maxiter )
        success=0; disp('Fail to converge within given number of iteration!');
        %subICmean_orth = O_inv * subICmean; %rotate final estimate of S
        return;
    end;
end;

% %rotate final estimate of S using inverse of orthonormalization matrix applied to A
% subICmean_orth = O_inv * subICmean;

end

