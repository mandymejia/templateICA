function [theta, subICmean, subICvar, success] = EM_easy (tempICmean, tempICvar, Y, theta0, C_matrix_diag, maxiter, epsilon)
% [theta, subICmean, subICvar, success] = EM_easy (tempICmean, tempICvar, Y, theta0, C_matrix_diag, maxiter, epsilon)
%
% EM algorithm for single-subject template ICA model (fast version with no nuisance/free ICs)
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
% theta0.A             : (QxQ) mixing matrix 
% theta0.nu0_sq        : (1x1) residual variance from first level
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

Q = size(Y, 1);
V = size(Y, 2);
itr = 1;
err= 1000; %set large initial value for difference between iterations
theta = theta0;
success = 1;
tempICvar(tempICvar < .00001) = .00001; %to prevent problems when inverting covariance

while err > epsilon
    [theta_new, subICmean, subICvar] = UpdateTheta_easy(Y, theta, C_matrix_diag, tempICmean, tempICvar);
    %subICmean, subICvar are mean and variance of subject-specific ICs
 
    A_vec = reshape(theta.A, [Q*Q,1]);
    A_vec_new = reshape(theta_new.A, [Q*Q,1]);
    err_A = norm(A_vec_new - A_vec)/norm(A_vec);
    err_nu = abs((theta_new.nu0_sq - theta.nu0_sq))/theta.nu0_sq;
    err = max(err_A, err_nu);

    fprintf('iteration %6.0f and the difference is  %6.6f for theta \n', itr, err);  
    clear A_vec A_vec_new err_A err_nu;
    theta = theta_new;
    itr = itr + 1;
    if ( itr > maxiter )
        success=0; 
        disp('Fail to converge within given number of iteration!');
        return;
    end;
end;

end

