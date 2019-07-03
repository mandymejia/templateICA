function [Probz_zlfix, idx] = Cal_Probzl (PostProbz, z_dict, q, m)
% Probz_zlfix = Cal_Probzl (PostProbz, z_dict, q, m)
% Calculating marginal prob P( z_q(v) = m | Y(v), theta )
% given m in {1,2,...,M}
% given q in {1,2,...,Q}
% i.e., the source of q th IC is fixed to m at voxel v
% M^(Q-1) possible combinations vector z(v) in total
idx = find(z_dict(q,:) == m);
Probz_zlfix = sum(PostProbz(idx)); 
end

