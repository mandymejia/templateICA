function [A, W, S, skew, iq] = icatb_calc_nuisanceICs(Y, Q, nvox)

% ICA is performed on the reduced data set
% Inputs:  1) Y: single-subject data matrix of size TxV
%          2) Q: number of ICs
%          3) nvox: number of voxels (V)
%
% Outputs: 1) A: mixing matrix
%          2) W: unmixing matrix
%          3) S: QxV matrix of component spatial maps
%          4) skew: tells you whether sign of component flipped so that most voxels positively contributing to it
%          5) iq: numerical vector of length Q; indicates quality of each component
%             based on repeated runs of Infomax algorithm (usually only use components with iq>.9)

if ~exist('statusHandle', 'var')
    statusHandle = [];
end

icatb_defaults;
global INDIVIDUAL_ICA_INDEX;
global GROUP_ICA_INDEX;
global WRITE_COMPLEX_IMAGES;
global NUM_RUNS_GICA;
global OPEN_DISPLAY_GUI;


% which_analysis = 1; run ICA once
which_analysis = 2; %run ICA multiple times with random initial conditions
icasso_opts.sel_mode = 'randinit';  % Options are 'randinit', 'bootstrap' and 'both'
icasso_opts.num_ica_runs = 5; % Number of times ICA will be run
% Most stable run estimate is based on these settings.
icasso_opts.min_cluster_size = 2; % Minimum cluster size
icasso_opts.max_cluster_size = icasso_opts.num_ica_runs; % Ma

icaAlgo = icatb_icaAlgorithm; % available ICA algorithms
% 1 means infomax, 2 means fastICA, etc.
algoVal = 1; % algorithm index

% selected ICA algorithm
algorithmName = deblank(icaAlgo(algoVal, :));

% Modality type
modalityType = 'fMRI';
dataTitle = 'Functional';
compSetFields = {'ic', 'tc'};

sesInfo.num_runs_gica = 1;


%% Open parallel mode
parallel_info.mode = 'serial';
parallel_info.num_workers = 5;

toolboxNames = ver;
parallelCluster = ~isempty(find(strcmpi(cellstr(char(toolboxNames.Name)), 'parallel computing toolbox') == 1));

if (strcmpi(parallel_info.mode, 'parallel'))
    statusHandle = [];
end

sesInfo.dataType = 'real';

%number of components to extract
numOfIC = Q;
% data = reshape(data,xdim*ydim*zdim,size(data,2));
mask_ind = 1:nvox;

ICA_Options = {};
ICA_Options = {'posact', 'off', 'sphering', 'on', 'bias', 'on', 'extended', 0};

% convert to cell
if isempty(ICA_Options)
    if ~iscell(ICA_Options)
        ICA_Options = {};
    end
end

% ICASSO

if (strcmpi(parallel_info.mode, 'serial') || parallelCluster)
    
    %%%%% Calculate PCA and Whitening matrix %%%%%
    % PCA
    [V, Lambda] = icatb_v_pca(Y, 1, numOfIC, 0, 'transpose', 'yes');
    
    % Whiten matrix
    [w, White, deWhite] = icatb_v_whiten(Y, V, Lambda, 'transpose');
    
    clear V Lambda;
    
    if strcmpi(parallel_info.mode, 'serial')
        sR = icatb_icassoEst(icasso_opts.sel_mode, Y, icasso_opts.num_ica_runs, 'numOfPC', numOfIC, 'algoIndex', algoVal, ...
            'dewhiteM', deWhite, 'whiteM', White, 'whitesig', w, 'icaOptions', ICA_Options);
    else
        % use PCT
        sR = icatb_parIcassoEst_cluster(icasso_opts.sel_mode, Y, icasso_opts.num_ica_runs, 'numOfPC', numOfIC, 'algoIndex', algoVal, ...
            'dewhiteM', deWhite, 'whiteM', White, 'whitesig', w, 'icaOptions', ICA_Options);
    end
    
    clear data w deWhite White;
    
    %%%% Visualization %%%%%%
    
    sR = icassoExp(sR);
    
    %%% Return icasso results  
    %IC quality
    iq = icassoResult(sR, numOfIC);
     
    minClusterSize = icasso_opts.min_cluster_size;
      
    maxClusterSize = icasso_opts.max_cluster_size;
    
    if (minClusterSize <= 1)
        minClusterSize = 2;
    end
    
    if (minClusterSize > icasso_opts.num_ica_runs)
        minClusterSize = icasso_opts.num_ica_runs;
    end
    
    [metric_Q, A, W, icasig] = getStableEstimates(sR, minClusterSize, maxClusterSize);
    
    clear sR;
    
    skew = zeros(1, numOfIC);
    
%     [A, W, icasig, skew] = changeSignOfComponents(A, icasig);
    S = icasig;
%    parallel_info.mode = 'serial';

end








