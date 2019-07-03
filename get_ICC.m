function [ICC, I2C2, var_W, var_U] = get_ICC(meas1, meas2, img_dim)
% [ICC, I2C2, var_W, var_U] = get_ICC(meas1, meas2, img_dim)
%
% INPUTS
% meas1 (n x p1 x ... x pq) : first set of measurements from n subjects on p = prod(pj) variables
% meas2 (n x p1 x ... x pq) : second set of measurements from n subjects on p = prod(pj) variables
% img_dim : which dimension AFTER THE FIRST corresponds to image locations (or other variables to collapse over for I2C2)
%			(e.g. if p1 is the number of voxels, then img_dim=1)
%			if img_dim=0, do not compute I2C2
%
% OUTPUTS
% ICC (p1 x ... x pq).  : intra-class correlation coefficient of each variable 
% var_W (p1 x ... x pq) : total variance (within-subject + between-subject) of each variable
% var_U (p1 x ... x pq) : within-subject variance of each variable

	dims = size(meas1);
	singlevardim = (numel(dims) == 2);
	n = dims(1);
	p = prod(dims(2:end));
	if(p > 1 & ~singlevardim) %if there's only one dimension representing variables, then the measurements are already vectors
        'reshaping data to matrix form'
        meas1_mat = reshape(meas1, [n, p]);
        meas2_mat = reshape(meas2, [n, p]);
	else
		meas1_mat = meas1;
		meas2_mat = meas2;
	end

	'computing total variance'
	var_W1 = var(meas1_mat);
	var_W2 = var(meas2_mat);
	var_W = (var_W1 + var_W2)/2; %total variance

	'computing within-subject variance'
	diffs_mat = meas1_mat - meas2_mat;
	var_U = (1/2)*var(diffs_mat); %within-subject variance

	'computing ICC'
	ICC = 1 - var_U./var_W; ICC(ICC < 0) = 0; ICC(ICC > 1) = 1; %ICC

	if(p > 1 & ~singlevardim) %if there's only one dimension representing variables, then the measurements were already vectors
		'reshaping data to array form'
		ICC = reshape(ICC, dims(2:end));
		var_W = reshape(var_W, dims(2:end));
		var_U = reshape(var_U, dims(2:end));
	end

	if(img_dim > 0)
		'computing I2C2'
		if(~singlevardim) I2C2 = 1 - sum(var_U, img_dim)./sum(var_W, img_dim); end
		if(singlevardim) I2C2 = 1 - sum(var_U)./sum(var_W); end
		I2C2(I2C2 < 0) = 0; I2C2(I2C2 > 1) = 1;
	else 
		I2C2 = NaN;
	end

end



