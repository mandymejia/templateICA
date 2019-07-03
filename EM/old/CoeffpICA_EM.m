function [theta, beta, z_mode, subICmean, subICvar, grpICmean, grpICvar, success] = CoeffpICA_EM (Y, X, theta0, C_matrix_diag, beta0, maxiter, epsilon1, epsilon2)
% [theta, beta] = CoeffpICA_EM (Y, X_mtx, theta0, beta0, N, T, q, p, m, V, maxiter, epsilon)
% Explicit form of EM algorith for the pICA model with interactions
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% beta        :  beta (k, l, v)  coefficients at voxel v,          p*q*V
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% theta: parameter space
% theta.A             : [A1,A2,...,AN]  T*q*N
% theta.sigma1_sq     : a scalar of sigma1_sq
% theta.sigma2_sq     : a vector of sigma21_sq,...,sigma2q_sq, 
%                              q variantions for each IC, q*1
% theta.miu3          : vector of miu_z, mq*1
% theta.sigma3_sq     : vector of sigma_z, mq*1
% theta.pi            : matrix of pi, mq*1
% beta0, theta0: initial guesses
X_mtx=X';
N=size(X,1);
p=size(X,2);
q=size(theta0.sigma2_sq, 1);
T=q;
m=2;
V=size(Y, 2);
itr = 1;
err1= 1000;
err2= 1000;
theta = theta0;
beta  = beta0;
success = 1;
while err1 > epsilon1 || err2 > epsilon2 
    %%[theta_new, beta_new] = UpdateThetaBeta (Y, X_mtx, theta, C_matrix_diag, beta, N, T, q, p, m, V);
    [theta_new, beta_new, z_mode, subICmean, subICvar, grpICmean, grpICvar, error] = UpdateThetaBeta (Y, X_mtx, theta, C_matrix_diag, beta, N, T, q, p, m, V);
    if(error == 1)
        success = 0;
        disp('Fail to converge due to calculate of p(z|y)!');
        return;
    end;
    [vec_theta_new, vec_beta_new] = VectThetaBeta ( theta_new, beta_new, p, q, V, T, N, m);
    [vec_theta, vec_beta]         = VectThetaBeta ( theta, beta, p, q, V, T, N, m);
    err1 = norm(vec_theta_new - vec_theta)/norm(vec_theta);
    err2 = norm(vec_beta_new  - vec_beta)/norm(vec_theta);  
    fprintf('iteration %6.0f and the difference is  %6.6f for theta and %6.6f for beta \n', itr, err1, err2);  
    %%%clear vec_theta_new vec_beta_new vec_theta vec_beta theta beta;
    theta = theta_new;
    beta  = beta_new;
    itr = itr + 1;
    if ( itr > maxiter )
        success=0; disp('Fail to converge within given number of iteration!');
        return;
    end;
end;
end

