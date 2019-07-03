function [H, D, U, sigma_sq, q] = dim_reduce(vectors, q)
	%USAGE
	%
	% [H, D, U, sigma_sq, q] = dim_reduce(vectors, q)
	% [H, D, U, sigma_sq, q] = dim_reduce(vectors, 0) to estimate q
	%
	% ARGUMENTS
	%
	% vectors	Data in row vectors (obs x variables)
	% q 		Number of dimensions to retain, or 0 if estimated

	p = size(vectors, 2); %number of variables
	n = size(vectors, 1); %number of observations

	% perform PCA
    [U_all, D_all] = pcamat(vectors, 1, n, 'off', 'off'); 
    [lambda,o] = sort(diag(D_all),'descend'); %sort eigenvalues

    %estimate dimensionality
    if(q==0)
    	[q, ~] = laplace_pca([], lambda, p, n);
    end

    %apply dimension reduction
    inds = o(1:q); % indices of largest q eigenvecs
    U = U_all(:,inds); % U matrix (n x q)
    D = diag(D_all(inds, inds)); % diagonals of D matrix (q x q)
    sigma_sq = mean(lambda((q+1):length(lambda))); % residual variance
    H = diag((D - sigma_sq).^(-1/2)) * U'; % prewhitening matrix

end
