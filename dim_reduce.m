function [Y_new, H, Hinv, D, U, sigma_sq, C_diag, Q] = dim_reduce(Y, Q)
	%USAGE
	%
	% [Y_new, H, Hinv, D, U, sigma_sq, C_diag] = dim_reduce(Y, Q)
	% [Y_new, H, Hinv, D, U, sigma_sq, C_diag, Q] = dim_reduce(Y, 0) to estimate Q
	%
	% ARGUMENTS
	%
	% Y (TxV) – BOLD data
	% Q 	  – Number of dimensions to retain, or 0 if to be estimated
    %
    % RETURNS
    %
    % Y_new (QxV) – dimension-reduced data 
    % H (QxT)     - dimension reduction / prewhitening matrix
    % Hinv (TxQ)  - reverse dimension reduction matrix
    % D (Qx1)     - diagonals of D matrix from SVD 
    % U (TxQ)     - U matrix from SVD
    % sigma_sq    - residual variance, based on last T-Q eigenvals
    % C_diag (Qx1)- diagonals of residual cov matrix in ICA model
    % Q           - dimensionality (specified or estimated)

	V = size(Y, 2); 
	T = size(Y, 1); 

	% perform PCA
    [U_all, D_all] = pcamat(Y, 1, T, 'off', 'off'); 
    [lambda,o] = sort(diag(D_all),'descend'); %sort eigenvalues

    %estimate dimensionality
    if(Q==0)
    	[Q, ~] = laplace_pca([], lambda, V, T);
    end

    %apply dimension reduction
    inds = o(1:Q); % indices of largest Q eigenvecs
    U = U_all(:,inds); % U matrix (T x Q)
    D = diag(D_all(inds, inds)); % diagonals of D matrix (Q x Q)
    sigma_sq = mean(lambda((Q+1):length(lambda))); % residual variance
    H = diag((D - sigma_sq).^(-1/2)) * U'; % prewhitening matrix (Q x T)
    Y_new = H * Y;
    
    %other things needed for ICA model
    C_diag = (D - sigma_sq).^(-1); % residual covariance
    Hinv = U * diag((D - sigma_sq)).^(1/2); %reverse dimension reduction matrix

end
