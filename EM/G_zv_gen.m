function G_zv = G_zv_gen(zv, M, Q) 
% G_zv = G_zv_gen(zv, M, Q) 
% Generates G_z matrix (Qx(MQ)) given vector z(v)
%
% INPUTS:
% zv - {zv1, zv2,..., zvQ}', zvi in {1,2,...,M} - one possible value of z(v)
% M - number of MoG components
% Q - number of ICs assumed to follow MoG distribution ('free' ICs)
%
% OUTPUT: 
% G_z - Qx(MQ) matrix given vector z(v)

	x=(1:Q)';
	y=(x-1)*M + zv;
	G_zv = zeros (Q, M*Q);
	G_zv( sub2ind(size(G_zv), x, y) ) = 1;

end
 