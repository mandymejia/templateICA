function [S_est_match, matched, corr_match] = match_components(S_est, S_true, sign_match)
	% USAGE
	%
	% [S, A, dat_ctr] = dual_reg(dat, S_grp)
    %
    % Matches estimated ICs to true ICs.  
    % Allows the number of estimated ICs to be mis-specified.
	%
	% ARGUMENTS
	%
	% S_est			Estimated independent components (VxQ)
	% S_true		True independent components (VxQ)
	% sign_match	1 - consider sign of correlation 
	%				0 - only consider magnitude of correlation (negative ok)
	%
	% OUTPUT
	%
	% S_est_match 	Re-ordered estimated independent components (VxQ)
	% corr_match	Correlation of true and estimated ICs after matching

    corr_mat = corr(S_true, S_est);
   	if(sign_match) 
		corr_mat(corr_mat < 0) = 0;
	else
		corr_mat = abs(corr_mat);
	end

    V = size(S_est, 1);
    Q = size(S_true, 2);
    if(V ~= size(S_true,1)) error('Different number of voxels in true and estimated S'); end

    matched = [];
    corr_match = [];
    while(max(max(corr_mat)) > 0) %while there are any more candidate matches

    	[I,J] = find(corr_mat == max(max(corr_mat)));
        corr_match = [corr_match, max(max(corr_mat))];
    	q_target = I; %index of the true IC with the strongest match to an estimated IC
    	q_match = J; %index of the estimated IC matched to the target
    	if(q_target ~= q_match)
    		S_match = S_est(:,q_match); %IC matched to target
    		S_switch = S_est(:,q_target); %the IC to switch with
    		S_est(:,q_target) = S_match;
    		S_est(:,q_match) = S_switch;
    	end
    	matched = [matched, q_target];
    	corr_mat = corr(S_true, S_est);

    	%eliminate already matched components from consideration
    	corr_mat(matched,:) = 0; 
    	corr_mat(:,matched) = 0;

    	%deal with negative correlations
    	if(sign_match) 
    		corr_mat(corr_mat < 0) = 0;
    	else
    		corr_mat = abs(corr_mat);
    	end

    end

    S_est_match = S_est;

end

