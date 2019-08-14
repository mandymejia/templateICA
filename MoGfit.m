%%%%% fit a Gaussin mixture
% Y is the data vector
% m is the number of components 
% Output includes mean, variance and pi
function [theta]=MoGfit(Y,m)
    Y = reshape(Y,prod(size(Y)),1);
    GMModel = gmdistribution.fit(Y,m,'Options',statset('MaxIter',1000));
    theta.miu(1:m, 1) = GMModel.mu;
    theta.sigma_sq(1:m, 1) = GMModel.Sigma;
    theta.pi(1:m, 1)  = GMModel.PComponents;
end

