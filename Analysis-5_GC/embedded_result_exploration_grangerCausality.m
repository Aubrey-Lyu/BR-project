%% SETUP
clear
% ------------ paths -----------------------
hhm_master_dir = '/media/dian/D/data/scripts/Github/HMM-MAR-master'; % adapt to yours

addpath(genpath(hhm_master_dir))
addpath(genpath('/home/dian/scripts/Github/Violinplot-Matlab-master'))
addpath(genpath('/media/dian/D/ProgramData/MATLAB_toolbox/hex_and_rgb_v1.1.1'))

dir_base = '/media/dian/D/data/Binocular_Rivalry/HMM_analyses';
cd(dir_base)
run_indx = 1;

%% what analyses


%% crucial parameters
%ROIs = {'rHP', 'lHP', 'lIPL', 'rIPL', 'PCC', 'PCU','ACC'};
conditions = {'BR_dominant';  'BR_mixed';
    'Rpl_smt_dominant';  'Rpl_smt_mixed';
    'Rpl_ins_dominant';  'Rpl_ins_mixed';
    }; % for epoch estimation, has to have same number of timepoints

Ncond = length(conditions);
cdt_code = eye(Ncond);

% scan files which you took data from; loop in the same order
load(fullfile(dir_base, 'data', 'CONDSCANS_for_voi_Y.mat'))

%% load data
%  key data
load(fullfile(dir_base,'data','data_DMN_V1-1000_100_TUDA.mat'));
load(fullfile(dir_base,'data','T_DMN-1000_100_TUDA.mat'));
% covariates
load(fullfile(dir_base,'data','Lat_before_DMN-1000_100_TUDA2.mat'));
load(fullfile(dir_base,'data','Lat_after_DMN-1000_100_TUDA2.mat'));
load(fullfile(dir_base,'data','Rsp_moments_DMN-1000_100_TUDA.mat'));
load(fullfile(dir_base,'data','CDT_code_DMN-1000_100_TUDA.mat'));
load(fullfile(dir_base,'data','previous_stable_dur_DMN-1000_100_TUDA.mat'));
% result data
%output_dir = [dir_base '/output_test/tuda_epoch_K' num2str(K) 'order' num2str(order) '_DMN-1000_100_meanact'];
output_dir = '/media/dian/D/data/Binocular_Rivalry/HMM_analyses/output_test/embedded_lag/DMN_V1_1000_100_K4_pca16';

cd(output_dir)

%% load result file
%load(fullfile(output_dir, ['run' num2str(run_indx) '.mat']));

%% key parameters
K = 4 ;% options.K;
Hz = 250; %options.Fs; order = options.order;

% if isfield(options, 'embeddedlags')
%     order = length(options.embeddedlags)-1;
% end

%% -------------  rewrite data  -----------------

Ts = T; % rewrite T
X  = [];
Y  = [];
T  = [];
Sub = [];
% prepare to rewrite data -->  cell to mat

disp('Loading data...')

Sn = length(Ts); Ts_new = cell(20,1);

for ss = 1:Sn
    ttrial = Hz*1; % just take one second %Ts{ss}(1); % assuming the ttrial is the same for all T
    xInd = []; n = 0; % counter for T{ss}.
    sub = ['S' sprintf('%02d%',ss)];
    
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
        %         if isempty(find(isnan(Lat_before{ss}(slice)))) && isempty(find(isnan(Lat_after{ss}(slice))))
        %             % rewrite data to have same length between trials
        xInd = [xInd, slice];
        n = n+1; % update counter
        Ts_new{ss}(n,1) = ttrial;
        %         end
        
        Sub = [Sub; sub];
    end
    % check inside variables
    disp(' ')
    disp(['Subj:   ' num2str(ss)])
    disp(['Ntrial: ' num2str(length(Ts{ss}))])
    disp(['N:      ' num2str(n)])
    
    data{ss} = data{ss}(xInd,:);
    X = [X; data{ss}];
    Y = [Y; Lat_before{ss}(xInd,:) Lat_after{ss}(xInd,:) previous_stable_dur{ss}(xInd,:) Rsp_moments{ss}(xInd,:) CDT_code{ss}(xInd,:)];
    T = [T; Ts_new{ss}];
    
end % for ss

disp('Data are loaded!')

%% --------------------------------------- %%

%%
for i = 1%:%Ncond
    disp(['Ncond = ' num2str(i)])
    condition = conditions{i};
    Ind = find( Y(:, 5+i)==1 );
    Y_  = Y(Ind, :);
    X_  = X(Ind, :);
    
    zeroMod_index = find(mod(Ind,ttrial) ==0);
    divisor       = floor(Ind/ttrial);
    iTT           = divisor(zeroMod_index);
    
    cond_sub = Sub(iTT,:);
    
    % gamma = []; vp = []; % vpath
    %     for i_t = 1:length(iTT)
    %         iT = iTT(i_t);
    %         gamma = [gamma; Gamma((sum(T(1:iT-1))-(iT-1)*order)+1 : (sum(T(1:iT-1))-(iT-1)*order)+ttrial-order, :)];
    %         %vp    = [vp;    vpath((sum(T(1:iT-1))-(iT-1)*order)+1 : (sum(T(1:iT-1))-(iT-1)*order)+ttrial-order, :)];
    %     end
    
    N   = size(Y_,1)/ttrial;
    x   = reshape(X_, [ttrial, N, 8]);
    Y0  = reshape(Y_, [ttrial, N, size(Y_,2)]);
    T_  = ones(N, 1) * ttrial;
    
end % for i = 1:Ncond

%% granger causality
addpath('/media/dian/D/data/scripts/matlab/granger_causality/granger_cause_1')
plot_dir = '/media/dian/D/data/Binocular_Rivalry/HMM_analyses/output_test/embedded_lag/DMN_V1_1000_100_K4_pca16/run1_analyses/fig';

max_x_lag = 15;
max_y_lag = 20;
alpha = 0.05; % significance level
H = 0; % if reject null hypothesis
ntrials = size(x,2);  %Define the number of trials.
firstYlag = 0;
b_len = max_x_lag+1+max_y_lag+(1-firstYlag);
use_best_x = 1;
use_best_y = 0;
% test ----------------------------
 pcu = squeeze(x(:,1,5));
v1  = squeeze(x(:,1,8));
% v1_cause_pcu
    [ F, c_v,  Fprob , ~,   dAIC, dBIC  , chosen_x_lag, chosen_y_lag, b ] ...
        =   granger_cause_1(pcu,v1, alpha, max_x_lag, 0, max_y_lag, 1, firstYlag ,  ...
        'PCu', 'V1' , 0 , plot_dir , ['trial ' num2str(itrial)],...
        ['Granger causality analysis for trial ' num2str(itrial)] , do_plot);
        
    v1_cause_pcu_temp = [F, Fprob, dAIC, dBIC, chosen_x_lag, chosen_y_lag, H, b'];
% -------------------------

pcu_cause_v1 = zeros(ntrials, length(v1_cause_pcu_temp));
v1_cause_pcu = zeros(ntrials, length(v1_cause_pcu_temp));


for itrial=1:ntrials  %Compute the spectra for each trial.
    pcu = normalize(squeeze(x(:,itrial,5)));
    v1  = normalize(squeeze(x(:,itrial,8)));
    
    if itrial <= 10; do_plot=0; else; do_plot=0; end
    % v1_cause_pcu
    [ F, c_v,  Fprob , ~,   dAIC, dBIC  , chosen_x_lag, chosen_y_lag, b ] ...
        =   granger_cause_1(pcu,v1, alpha, max_x_lag, use_best_x, max_y_lag, use_best_y, firstYlag ,  ...
        'PCu', 'V1' , 0 , plot_dir , ['trial ' num2str(itrial)],...
        ['Granger causality analysis for trial ' num2str(itrial)] , do_plot);
    
    if F>c_v; H=1; end
    if length(b) < b_len
        display(['Trial ' num2str(itrial) ' for V1-PCu only have '  num2str(length(b)) ' betas.'])
        bbb = ones(b_len,1)*NaN; 
        bbb(1:chosen_x_lag+1) = b(1:chosen_x_lag+1);
        bbb(max_x_lag+2 : end) = b((chosen_x_lag+2): (chosen_x_lag+1+max_y_lag+(1-firstYlag)));
        b = bbb;
    end
        
    v1_cause_pcu(itrial,:) = [F, Fprob, dAIC, dBIC, chosen_x_lag, chosen_y_lag, H, b'];
    
    clear F c_v  Fprob  Fprob_Corrected  dAIC dBIC   chosen_x_lag chosen_y_lag b
    
    % pcu_cause_v1
    [ F, c_v,  Fprob , Fprob_Corrected,   dAIC, dBIC  , chosen_x_lag, chosen_y_lag, b ] ...
        =   granger_cause_1(v1, pcu, alpha, max_x_lag, use_best_x, max_y_lag, use_best_y, firstYlag ,  ...
        'V1', 'PCu', 0 , plot_dir , ['trial ' num2str(itrial)], ...
        ['Granger causality analysis for trial ' num2str(itrial)] , do_plot);
    
    if F>c_v; H=1; end
    if length(b) < b_len
        display(['Trial ' num2str(itrial) ' for PCu-V1 only have '  num2str(length(b)) ' betas.'])
        bbb(1:chosen_x_lag+1) = b(1:chosen_x_lag+1);
        bbb(max_x_lag+2 : end) = b((chosen_x_lag+2): (chosen_x_lag+1+max_y_lag+(1-firstYlag)));
        b = bbb;
    end
    pcu_cause_v1(itrial,:) = [F, Fprob, dAIC, dBIC, chosen_x_lag, chosen_y_lag, H, b'];
    
end
save('/media/dian/D/data/Binocular_Rivalry/HMM_analyses/output_test/embedded_lag/DMN_V1_1000_100_K4_pca16/run1_analyses/v1_cause_pcu.mat', 'v1_cause_pcu')
save('/media/dian/D/data/Binocular_Rivalry/HMM_analyses/output_test/embedded_lag/DMN_V1_1000_100_K4_pca16/run1_analyses/pcu_cause_v1.mat', 'pcu_cause_v1')
% 
% %% Filter trials based on state life times
% filter_k4lifetime = zeros(length(lifetimes),1);
% k4_life = zeros(length(lifetimes),1);
% for i = 1:length(lifetimes)
%     life = sum(lifetimes{i, 4});
%     k4_life(i) = life;
%     if life > 275*0.6
%         filter_k4lifetime(i) = 1;
%     end
% end
% 
% %% behav relationship
% dat = pcu_cause_v1;
% mask = find(filter_k4lifetime==1 & dat(:,4)<=0);
% 
% xx = log(-1*dat(mask,4));
% %xx = dat(mask,1);
% yy = Y0(1,mask,3);
% xx = xx(~isnan(yy));
% yy = yy(~isnan(yy));
% close all; figure; scatter(xx,yy)
% mdl = fitglm(xx,yy);
% 
% %% bhav relationship -- within subj
% 
% % get subject glm fitting
% B = [];
% dat = pcu_cause_v1;
% 
% close all
% for ss = 1:20
%     sub = ['S' sprintf('%02d', ss)];
%     mask = find(strcmp(cellstr(cond_sub), sub) & ...
%         dat(:,4)<-50 & filter_k4lifetime==1);
%     vy = log(Y0(1,mask,3));
%     vx = dat(mask,4);
%     vy = vy(~isnan(vy));
%     vx = vx(~isnan(vy));
%     if length(vx) >10
%         b = glmfit(vx, vy);
%         yfit = b(1) + vx*b(2);
%         figure; scatter(vx, vy); hold on; plot(vx, yfit,'--'); hold off
%         B = [B; b(2)];
%     else
%         continue
%     end
% end
% 
% [h,p]=ttest(B);



