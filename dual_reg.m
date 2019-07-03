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

	%center timeseries data across space and time and standardize scale
	dat = dat - mean(dat); %center each voxel time series (centering across time)
	dat_t = dat'; %transpose image matrix
  	sig = sqrt(mean(var(dat_t))); %variance across image, averaged across time, square root to get SD
  	dat = (dat_t - mean(dat_t))'; %center each image (centering across space)
  	dat = dat./sig; %standardize by global SD
  	dat_ctr = dat;

  	%center group ICs over voxels
	S_grp = (S_grp' - mean(S_grp'));
	S_grp_t = S_grp';

	%estimate A (IC timeseries)
	A = dat_ctr * S_grp * inv(S_grp_t * S_grp);  
	%fix scale of timeseries
	sd_A = std(A); 
	D = diag(1./sd_A);
  	A = A * D;

	%estimate S (IC maps)
	S = inv(A' * A) * (A' * dat_ctr);

end
