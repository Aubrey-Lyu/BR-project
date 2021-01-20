%% analysis: HMM validation on different runs, all channel, train HMM

% The version two will take the pre-onset time window: [-1 0.1]
% Note: data have been reorganised in the HMM_meta_output.mat file
% for example, the 'data_DMN_V1-1000_100_TUDA.mat' used in this original
% analysis was called "X" in the uploaded HMM_meta_output.mat file

tic
%% what analyses?
organise_data = 1;
run_HMM       = 1;
use_stochastic =1;
loopin = randperm(20); % permute subjs order

% =========  run on half of the subjects ===================
for half = 1:2

%% define variable
ROIs = {'rHP', 'lHP', 'lIPL', 'rIPL', 'PCU', 'PCC', 'ACC', 'V1'};
HMM_dir = '/data/dl577/BR/HMM_analyses';
dir_base = '/data/dl577/BR/HMM_analyses';
cd(dir_base)

%% organise data
if organise_data
     % load key data
    load(fullfile(dir_base,'data','data_DMN_V1-1000_100_TUDA.mat'));
    load(fullfile(dir_base,'data','T_DMN-1000_100_TUDA.mat'));
    % load covariates
    load(fullfile(dir_base,'data','Lat_before_DMN-1000_100_TUDA.mat'));
    load(fullfile(dir_base,'data','Lat_after_DMN-1000_100_TUDA.mat'));
    load(fullfile(dir_base,'data','Rsp_moments_DMN-1000_100_TUDA.mat'));
    load(fullfile(dir_base,'data','CDT_code_DMN-1000_100_TUDA.mat'));
    
    % rewrite data
    Ts = T; % rewrite T
    X  = [];
    Y  = [];
    T  = [];
    
    % prepare to rewrite data -->  cell to mat
    
    disp('Loading data...')
    
    Sn = length(Ts); Ts_new = cell(20,1);

    for ss = loopin((1:Sn/2) + (Sn/2)*(half-1))
        ttrial = Ts{ss}(1); % assuming the ttrial is the same for all T
        xInd = []; n = 0; % counter for T{ss}.
        
        for i_ntrial = 1:length(Ts{ss})
            
            % get rid of NaN value; this problem only exits in TUDA analyses
            % determine if in this trial, Lat is NaN, if it is, skip
            % this trial
            if i_ntrial == 1
                slice = 1:ttrial;
            else
                slice = sum(Ts{ss}(1:i_ntrial-1))+1 : sum(Ts{ss}(1:i_ntrial-1)) + ttrial;
            end
            %--------------------------------------------------------------------------------------
           % if isempty(find(isnan(Lat_before{ss}(slice)))) && isempty(find(isnan(Lat_after{ss}(slice))))
                % rewrite data to have same length between trials
                xInd = [xInd, slice];
                n = n+1; % update counter
                Ts_new{ss}(n,1) = ttrial;
           % end
        end
        % check inside variables
        disp(' ')
        disp(['Subj:   ' num2str(ss)])
        disp(['Ntrial: ' num2str(length(Ts{ss}))])
        disp(['N:      ' num2str(n)])
        
        data{ss} = data{ss}(xInd,:);
        X = [X; data{ss}];
        Y = [Y; Lat_before{ss}(xInd,:) Lat_after{ss}(xInd,:) Rsp_moments{ss}(xInd,:) CDT_code{ss}(xInd,:)];
        T = [T; Ts_new{ss}];
        
    end % for ss
    
    disp('Data are loaded!')
    
end

%% ========     run HMM: embedded lagging   =======================
%{
The idea is that we model the main principal components (eigenvectors) of
a Gaussian distribution defined not only over space but also over a window
of time around the point of interest. This approach, inspired by and related
to the theory of Gaussian processes, can capture both spatial and spectral
 properties without overfitting. The limitation is that the use of dimensionality
reduction over this (space x time) by (space x time) matrices that represent
the states makes the approach to be relatively blind to the high frequencies.
%}
% =================================================================
if run_HMM
    hhm_master_dir = '/data/dl577/scripts/matlab/GitHub/HMM-MAR-master'; % adapt to yours
    addpath(genpath(hhm_master_dir))
    disp('Preparing to run HMM - embedded lag model...')
    K = 4; % no. states
    order = 0;
    Hz    = 250;
    N     = length(T); % num of subj/sess
    disp(['K=' num2str(K) ', Hz=' num2str(Hz)])
    % define option
    options = struct();
    options.K = K; % number of states
    options.order = 0; % no autoregressive components
    options.embeddedlags=-7:7;
    options.covtype='full';
    options.zeromean = 1;
    options.pca = 16; % reduce dimension
    options.verbose = 1;
    options.standardise = 1;
    options.inittype = 'HMM-MAR';
    options.cyc = 500;
    options.initcyc = 10;
    options.initrep = 3;
    options.Fs = Hz;
    % stochastic options
    if use_stochastic
        options.BIGNbatch = 8;
        options.BIGtol = 1e-7;
        options.BIGcyc = 1500;
        options.BIGundertol_tostop = 5;
        options.BIGforgetrate = 0.7;
        options.BIGbase_weights = 0.9;
    end
    %% run HMM
    [hmm, Gamma, ~, vpath] = hmmmar(X,T,options);
    output_dir = fullfile(HMM_dir, 'embedded_lag', ['DMN_V1_1000_100_K' num2str(K) '_pca16']);

    mkdir(output_dir)
    
    result_file = dir(fullfile(output_dir, ['run*_half' num2str(half) '.mat']));
    if isempty(result_file)
        run_indx = 0;
    else
        result_name = result_file(end).name;
        run_indx = strsplit(result_name, '_');
        run_indx = strsplit(run_indx{1}, 'run');
        run_indx = str2double(run_indx{2});
    end
   
    save(fullfile(output_dir, ['run' num2str(run_indx+1) '_half' num2str(half) '.mat']), 'Gamma','vpath','hmm', 'options')
    disp(['HMM done for half' num2str(half) '.'])
    disp('==================================')
    disp(' ')
end
end % half
toc