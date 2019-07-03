function [S, A, dat_ctr] = dual_reg(dat, S_grp)
	% USAGE
	%
	% [S, A, dat_ctr] = dual_reg(dat, S_grp)
	%
	% ARGUMENTS
	%
	% dat		Subject-level fMRI data (TxV) 
	% S_grp 	Group-level independent components (QxV)
	%
	% OUTPUT
	% S 		Subject-level independent components (QxV)
	% A 		Subject-level mixing matrix (TxQ)
	% dat_ctr	The row- and column- centered fMRI data

	V = size(dat, 2); %number of voxels
	T = size(dat, 1); %number of time points
	Q = size(S_grp, 1); %number of ICs

	if(T > V) warning('More time points than voxels. Are you sure?'); end
	if(Q > V) warning('More ICs than voxels. Are you sure?'); end
	if(Q > T) warning('More ICs than time points. Are you sure?'); end
	if(V ~= size(S_grp, 2)) error('The number of voxels in dat and S_grp must match'); end

	dat_ctr = zscore(dat); %center and scale each voxel time series
	dat_ctr = (dat_ctr' - mean(dat_ctr'))'; %center each time point 
	S_grp = (S_grp' - mean(S_grp'));
	S_grpt = S_grp';
	A = dat_ctr * S_grp * inv(S_grpt * S_grp);
	S = inv(A' * A) * (A' * dat_ctr);

end
