function [theta_new, z_mode, subICmean, subICvar, error, zero_post_z] = UpdateTheta_subspace(Y, theta, C_matrix_diag, tempICmean, tempICvar)
% [theta_new, z_mode, subICmean, subICvar, error, zero_post_z] = UpdateTheta_subspace(Y, theta, C_matrix_diag, tempICmean, tempICvar)
%
% %%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%
%
% Y (QxV matrix)    : PCA-prewhitened fMRI data 
% theta (structure) : current parameter estimates
% theta.A           : (QxQ) mixing matrix 
% theta.nu0_sq      : (1x1) residual variance from first level
% theta.miu         : (M(Q-L)x1) miu_z (MoG means)
% theta.sigma_sq    : (M(Q-L)x1) sigma^2 (MoG variances)
% theta.pi          : (M(Q-L)x1) pi (MoG component probabilities)
% C_matrix_diag (Qx1) : diagonal elements of matrix proportional to residual variance.  
% tempICmean (LxV matrix) : mean of each IC in template
% tempICvar (LxV matrix)  : between-subject variance of each IC in template
%
% %%%%%%%%%%%%%%% OUTPUTS %%%%%%%%%%%%%%%
%
% theta_new (structure) - updated parameter estimates
% z_mode - most likely MoG component membership for each voxel & IC
% subICmean - estimates of subject-level ICs
% subICvar - variance of subject-level ICs (for inference)
% error - flag indicating convergence (0) or not (1)
% zero_post_z - number of voxels for which posterior probability of z is zero within subspace

% First calculate the conditional probability of ICs given the data and latent sorce state

error = 0;  %%no error

V = size(Y, 2); %number of voxels
Q = size(theta.A, 1); %number of total ICs
L = size(tempICmean,1); %number of template ICs
M = size(theta.miu, 1)/(Q-L); %number of MoG components
Q2 = Q-L; %number of free ICs

theta_new.A         = zeros(Q, Q);
theta_new.nu0_sq    = 0;
theta_new.miu       = zeros(M,Q2);   %pi, miu, sigma in the order of miu_l1,...,miu_lm, l=1:q
theta_new.sigma_sq  = zeros(M,Q2);
theta_new.pi        = zeros(M,Q2);

A_ProdPart1         = zeros(Q, Q);  %first part of the product format for Ai in M-step
A_ProdPart2         = zeros(Q, Q);  %second part of the product format for Ai in M-step

z_mode = zeros(V, 1);

subICmean = zeros(Q, V); %subject-specific ICs mean
subICvar = zeros(Q, Q, V); %subject-specific ICs var

A   = theta.A ; 
A1 = A(:,1:L);
A2 = A(:,((L+1):Q));
C_inv = diag(1./C_matrix_diag); %QxQ
Sigma0     =   diag(C_matrix_diag.*theta.nu0_sq); %first-level covariance matrix
Sigma0_inv = diag(1./diag(Sigma0));  %first-level inverse covariance matrix

%Ying's Approach: Background MoG component = 1
%Assume at most one IC is activated (MoG != 1)

% %%%%% dictionary for the z(v)s: 

% % 1 2 ... M 1 ... 1 ... 1 ... 1
% % 1 1 ... 1 2 ... M ... 1 ... 1
% % . . ... . 1 ... 1 ... 1 ... 1
% % . . ... . . ... . ... . ... .
% % . . ... . . ... . ... . ... .
% % 1 1 ... 1 1 ... 1 ... 2 ... M

% % Each column is a possible value of z(v), 1 + (M-1)*Q2 columns

% %create z_dict matrix
% firstcol = ones(Q2,1); %first column
% %remaining columns
% torep = [2:M, ones(1,(M-1)*Q2)]; %first line, plus first block of second line
% toreshape = [repmat(torep, [1,Q2-1]), 2:M ]; %entire matrix except for first column, strung into a vector of size Q2*Q2*(M-1)
% othercols = reshape(toreshape, [(M-1)*Q2,Q2])'; 
% z_dict = [firstcol,othercols];

% My Approach: signal MoG component = M 
% Assume at most one IC is activated (MoG == M)
% Should be equivalent to Ying's if M=2 
% Otherwise, not as computationally advantageous

% Let the Mth MoG component represent "signal", first 1:(M-1) are background

% NO IC WITH z_q(v)=M
z_dict = zeros(Q2,(M-1)^Q2); 
for zi = 1:(M-1)^Q2
    z_dict(:,zi) = z_gen(zi-1, M-1, Q2);
end

% EXACTLY ONE IC WITH z_q(v)=M

%enumerate all background possibilities for (Q2-1) IC's 
num_bgd = (M-1)^(Q2-1);
z_dict_bgd = zeros(Q2-1,num_bgd);
for zi = 1:num_bgd
    z_dict_bgd(:,zi) = z_gen(zi-1, M-1, Q2-1);
end

%set qth IC to M, other (Q2-1) IC's background 
for q = 1:Q2
    z_dict_q = zeros(Q2, num_bgd);
    z_dict_q(q,:) = M; %set z_q(v)=M
    notq = setdiff(1:Q2,q);
    z_dict_q(notq,:) = z_dict_bgd; %set z_q'(v) to all background possibilities
    z_dict = [z_dict, z_dict_q];
end

% CREATE G_zv_dict ARRAY
size_dict = size(z_dict,2);
G_zv_dict = zeros(Q2, M*Q2, size_dict);
for zi = 1:size_dict
    G_zv_dict(:,:,zi) = G_zv_gen(z_dict(:,zi), M, Q2); % Q2 x (M*Q2) for each value of z(v)
end



%first part of MVN mean for computation of Probzv

%Quantities needed for posterior probabilities of z(v) & posterior moments of s(v)|z(v)
miu_s1v_condz_part2a = (A1' * C_inv * Y)./theta.nu0_sq; %mu = Cov*(part2a + part2b)
miu_s2v_condz_part2a = (A2' * C_inv * Y)./theta.nu0_sq; %mu = Cov*(part2a + part2b)
cov_sv_condz_part1 = (A' * C_inv * A)./theta.nu0_sq;
mvn_mean_part1 = A1 * tempICmean; %(QxL)*(LxV)

%LOOPING OVER ALL VOXELS (COULD USE A SAMPLE TO START)
zero_post_z = 0;
for v = 1:V

    %strcat(num2str(v),'/',num2str(V))

    %second-level variance of first L components
    Sigma_v = diag(tempICvar(:,v));
    Sigma_v_inv = diag(1./tempICvar(:,v));

    %E-Step: generating all moments for updating 
    Y_v  = Y(:, v);        

    %%%%% Save posterior probabilities of z(v) & posterior moments of s(v)|z(v)
    logProbzv = zeros(size_dict, 1);
    miu_sv_condz     = zeros(Q, 1, size_dict);
    miu_sv_svT_condz = zeros(Q, Q, size_dict);

    %%%%% Precompute quantities needed for posterior moments of s(v)|z(v)
    miu_s1v_condz_part2 = miu_s1v_condz_part2a(:,v) + Sigma_v_inv * tempICmean(:,v);
    cov_sv_condz_part2 = blkdiag(Sigma_v_inv, zeros(Q2));
    mvn_cov_part1v = A1 * Sigma_v * A1' + Sigma0;

    %LOOPING OVER THE PROBABILITY SPACE OF Z (WITHIN SUBSPACE) (2 seconds for M=3, Q2=10)
    logp_z = zeros(size_dict,1); %log prior probability
    for zi = 1:size_dict 

        % Grab mu_zv, sigma^2_zv AND pi_zv
        G_zv = G_zv_dict(:,:,zi); %multiplier matrix to get miu_zv, sigma^2_zv * pi_zv for current value of zv
        miu_zv = G_zv * theta.miu;
        sigma_sq_zv = G_zv * theta.sigma_sq;
        pi_zv      = G_zv * theta.pi;
        Sigma_zv = diag(sigma_sq_zv); %D_{z(v)}
        Sigma_zv_inv = diag(1./sigma_sq_zv);

        % Compute MVN density at y_v           
        mvn_mean_zv = mvn_mean_part1(:,v) + A2 * miu_zv;
        mvn_cov_zv = mvn_cov_part1v + (A2 * Sigma_zv * A2');

        % Compute P(z(v)|y(v),theta) on log scale
        % P(z(v)|y(v),theta) \propto prod(pi_z) * mvnpdf(y_v',mvn_mean, mvn_cov) ; 
        logp_z(zi) = sum(log(pi_zv)); %prod(pi_zv)
        logp_zv_mvn  =  log(mvnpdf(Y_v', mvn_mean_zv', mvn_cov_zv));
        logProbzv(zi) = logp_z(zi) + logp_zv_mvn; %correct up to a constant, compute and divide by denominator after loop

        % Compute posterior moments of s(v)|z(v)
        % covariance
        cov_sv_condz_part2_zv = cov_sv_condz_part2;
        cov_sv_condz_part2_zv((L+1):Q,(L+1):Q) = Sigma_zv_inv; %second block different for each value of z(v)
        cov_sv_condz = inv(cov_sv_condz_part1 + cov_sv_condz_part2_zv);
        % first moment
        miu_s2v_condz_part2 = miu_s2v_condz_part2a(:,v) + Sigma_zv_inv * miu_zv;
        miu_sv_condz(:,:,zi) = cov_sv_condz * [miu_s1v_condz_part2; miu_s2v_condz_part2];
        % second moment
        miu_sv_svT_condz(:,:,zi) = miu_sv_condz(:,:,zi) * miu_sv_condz(:,:,zi)' + cov_sv_condz;
    end

    %miu_s2v_condz = miu_sv_condz((L+1):Q,:,:);
    %miu_s2v_sv2T_condz = miu_sv_svT_condz((L+1):Q,(L+1):Q,:);

    % SCALE P(Z|Y,THETA) TO SUM TO 1

    if(max(logProbzv) > -Inf) %posterior probabilities not all zero
        %multiply probabilities by a constant to avoid numerical issues
        %(add constant on the log scale)
        logK = -1*min(logProbzv(logProbzv > -Inf));
        Probzv = exp(logProbzv + logK);
        PostProbz = Probzv./sum(Probzv);    
    else %use prior instead of posterior
        prior_z = exp(logp_z);
        PostProbz = prior_z./sum(prior_z);
        zero_post_z = zero_post_z + 1;
    end
    
    %z_mode(v) = find(PostProbz == max(PostProbz));
    
    if( sum(isnan(PostProbz)) > 0)    %%%%% Handling errors when evaluating mvnpdf
        disp('Error in calculating Pr(z(v)|Y(v),theta) !!!');
        error = 1;
        return;
        %flag = flag + 1;
        %vnum(flag) = v;
    end;

    % COMPUTE E[s2|y,theta] and E[s2*s2'|y,theta] (matrix algebra is ~5x faster than Cal_SvMmt function)
    miu_sv = squeeze(miu_sv_condz) * PostProbz;
    miu_sv_svT = reshape(reshape(miu_sv_svT_condz, [Q*Q, size_dict]) * PostProbz, [Q,Q]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% E STEP COMPLETE, BEGIN M STEP %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     subICmean(:,v) = miu_sv;
     subICvar(:,:,v) = miu_sv_svT - miu_sv * miu_sv';
     
     %%% UPDATE A
     %sum part1 and part2 over voxels, then multiply at the end to update A
     A_ProdPart1 = A_ProdPart1 + Y_v * miu_sv';
     A_ProdPart2 = A_ProdPart2 + miu_sv_svT;
     
     %%% UPDATE nu_0^2
     theta_new.nu0_sq = theta_new.nu0_sq + (Y_v' * C_inv * Y_v) - 2*(Y_v' * C_inv * A * miu_sv) + trace(A' * C_inv * A * miu_sv_svT);

     % UPDATE MOG PARAMETERS (iteratively sum over voxels)
     for q2 = 1:Q2

        miu_svq_condz = squeeze(miu_sv_condz((L+q2),:,:)); 
        miu_svq2_condz = squeeze(miu_sv_svT_condz((L+q2),(L+q2),:));

         for m = 1:M
               [prob_zqm, idx] = Cal_Probzl (PostProbz, z_dict, q2, m); %marginalize over q'\ne q to get Pr(z_q(v)=m|y,theta), idx indexes the elements of the dictionary that have z_q==m
               [miu_svq_fixzq, miu_svq2_fixzq] = ...
                               Cal_SvMmtCondzl_subspace( miu_svq_condz , miu_svq2_condz, PostProbz, idx, M, Q2, q2);
               theta_new.pi(m,q2)       =  theta_new.pi(m,q2) + prob_zqm;    %%%%% V*pi in fact
               theta_new.miu(m,q2)      =  theta_new.miu(m,q2) + prob_zqm * miu_svq_fixzq;
               theta_new.sigma_sq(m,q2) =  theta_new.sigma_sq(m,q2) + prob_zqm * miu_svq2_fixzq;
         end;
     end; 
	 		 
     %clear PostProbz miu_sv_all_condz miu_sv_svT_all_condz ;

end; %FINISH LOOP OVER VOXELS

clear miu_sv miu_sv_svT miu_sv1 miu_sv2 miu_s1v_s1vT

%%% UPDATE PARAMETERS
theta_new.nu0_sq = 1/(V*Q) * theta_new.nu0_sq;
theta_new.A = A_ProdPart1 * inv(A_ProdPart2); %new Ai
theta_new.A = theta_new.A * real(inv((theta_new.A' * theta_new.A)^(1/2))); % symmetric orthogonalization
% sd_A = std(Hinv * theta_new.A); %determine scale of A after reverse-prewhitening
% theta_new.A = theta_new.A * diag(1./sd_A); %standardize scale

%%% UPDATE MOG PARAMETERS

theta_new.pi        = 1/V * theta_new.pi;                % new pi 
theta_new.miu      = (1/V) * theta_new.miu./theta_new.pi;                         % new miu
theta_new.sigma_sq = (1/V) * theta_new.sigma_sq./theta_new.pi - theta_new.miu.^2;     % new sigma^2

theta_new.pi = reshape(theta_new.pi, [M*Q2,1]);
theta_new.miu = reshape(theta_new.miu, [M*Q2,1]);
theta_new.sigma_sq = reshape(theta_new.sigma_sq, [M*Q2,1]);






