function [k,p] = laplace_pca(data, e, d, n, alpha, beta)
% LAPLACE_PCA   Estimate latent dimensionality by Laplace approximation.
%
% k = LAPLACE_PCA([],e,d,n) returns an estimate of the latent dimensionality
% of a dataset with eigenvalues e, original dimensionality d, and size n.
% LAPLACE_PCA(data) computes (e,d,n) from the matrix data 
% (data points are rows)
% [k,p] = LAPLACE_PCA(...) also returns the log-probability of each 
% dimensionality, starting at 1.  k is the argmax of p.

% Written by Tom Minka; based on method described in
%     Minka, T. (2000).
%     Automatic choice of dimensionality for PCA.
%     Technical Report 514, MIT. 
% MBN updated to be compatible with Matlab2016b

if ~isempty(data)
  [n,d] = size(data);
  m = mean(data);
  data0 = data - repmat(m, n, 1);
  e = svd(data0,0).^2/n;
end
if nargin < 5
  alpha = 0;
  beta = 0;
end
e = e(:);
% break off the eigenvalues which are identically zero
ind = find(e < eps);
e(ind) = [];

% logediff(i) = sum_{j>i} log(e(i) - e(j))
logediff = zeros(1,length(e));
for i = 1:(length(e)-1)
  j = (i+1):length(e);
  logediff(i) = sum(log(e(i) - e(j))) + (d-length(e))*log(e(i)); %where does the 2nd part come from in eqn (78)
end
cumsum_logediff = cumsum(logediff);

n1 = n-1+alpha;
ehat = (e + beta/n)*n/n1;
inve = 1./ehat;
% invediff(i,j) = log(inve(i) - inve(j))  (if i > j)
%               = 0                       (if i <= j)
invediff = repmat(inve,1,length(e)) - repmat(inve',length(e),1);
invediff(invediff <= 0) = 1;
invediff = log(invediff);
% cumsum_invediff(i,j) = sum_{t=(j+1):i} log(inve(t) - inve(j))
cumsum_invediff = cumsum(invediff,1);
% row_invediff(i) = sum_{j=1:(i-1)} sum_{t=(j+1):i} log(inve(t) - inve(j))
% row_invediff = row_sum(cumsum_invediff);
row_invediff = sum(cumsum_invediff, 2);
% row_invediff(k) = sum_{i=1:(k-1)} sum_{j=(i+1):k} log(inve(j) - inve(i))

loge = log(ehat);
cumsum_loge = cumsum(loge);

cumsum_e = cumsum(ehat);

dn = length(e);
kmax = length(e)-1;
%dn = d;
%kmax = min([kmax 15]);
ks = 1:kmax;
% the normalizing constant for the prior (from James)
% sum(z(1:k)) is -log(p(U))
z = log(2) + (d-ks+1)/2*log(pi) - gammaln((d-ks+1)/2);
cumsum_z = cumsum(z);
for i = 1:length(ks)
  k = ks(i);
  %e1 = e(1:k);
  %e2 = e((k+1):length(e));
  %v = sum(e2)/(d-k);
  %v = (cumsum_e(end) - cumsum_e(k))/(length(cumsum_e)-k);
  v = (cumsum_e(end) - cumsum_e(k))/(d-k);
  p(i) = -n1/2*cumsum_loge(k) + (-n1*(d-k)/2)*log(v);
  p(i) = p(i) - cumsum_z(k) - k/2*log(n1);
  % compute h = logdet(A_Z)
  h = row_invediff(k) + cumsum_logediff(k); %Az_part_1 = cumsum_logediff(k), Az_part_2 includes row_invediff(k)
  % lambda_hat(i)=1/v for i>k
  h = h + (d-k)*sum(log(1/v - inve(1:k))); %Az_part_2 includes (d-k)*sum(log(1/v - inve(1:k)))
  m = d*k-k*(k+1)/2;
  h = h + m*log(n); %Az_part_3 is m*log(n)
  p(i) = p(i) + (m+k)/2*log(2*pi) - h/2;
  % missing terms added August 21 2008
  p(i) = p(i) + 1.5*k*log(2);
  p(i) = p(i) - 0.5*log(d-k);
  if alpha > 0
    ck = alpha*(d-k)/2*log(beta*(d-k)/2) - gammaln(alpha*(d-k)/2) + ...
	k*(alpha/2*log(beta/2) - gammaln(alpha/2));
    p(i) = p(i) - n1*d/2 - 0.5*log(n1) + ck;
  end
end
[pmax,i] = max(p);
k = ks(i);
  
v0 = cumsum_e(end)/length(cumsum_e);
p0 = -n1*d/2*log(v0) - 0.5*log(d);
%p0 = -n1*d/2*log(v0)
%p0 = -alpha*d/2*log(beta/alpha);
%p0 = -alpha*d/2*log(beta*d/2) + alpha*d/2*log(alpha*d/2);
if alpha > 0
  % must inline ck and put at end
  p0 = p0 - n1*d/2 - 0.5*log(n1) + alpha*d/2*log(beta*d/2) - gammaln(alpha*d/2);
end
if p0 >= pmax
  k = 0;
end
p = [p0 p];
