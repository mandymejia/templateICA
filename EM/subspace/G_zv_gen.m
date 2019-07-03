function G_zv = G_zv_gen(zv, m, q) 
% G_zv = G_zv_gen(zv, m, q) 
% zv={zv1, zv2,..., zvq}', zvi in {1,2,...,m}
% m: # of independent source of Gaussian Dist'n
% q: # of ICs
% G_z matrix generator given vector z(v): q*mq
 x=(1:q)';
 y=(x-1)*m + zv;
 G_zv = zeros (q, m*q);
 G_zv( sub2ind(size(G_zv), x, y) ) = 1;
end
 