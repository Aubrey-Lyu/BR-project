% % hhm_master_dir = '/data/dl577/scripts/matlab/GitHub/HMM-MAR-master'; % adapt to yours
% % addpath(genpath(hhm_master_dir))
% % 
% % dir_base   = '/data/dl577/BR/HMM_analyses';
% % result_dir = '/data/dl577/BR/HMM_analyses/embedded_lag/DMN_V1_1000_100_K4_pca16';
% % ROIs = {'rHP', 'lHP', 'lIPL', 'rIPL', 'PCU','PCC', 'ACC', 'V1'};
% % 
% % load(fullfile(result_dir, 'run1.mat'));
% % outputfile = [result_dir '/run1_analyses/hmm_analysis.mat'];
% % mkdir(fullfile(result_dir, 'run1_analyses'))
% % cd(fullfile(result_dir, 'run1_analyses'))
% % mkdir('fig')
% % 
% % %% get auto covariance matrix
% % close all
% % for i = 1:options.K
% %     figure; imagesc(getAutoCovMat(hmm,i)); colorbar
% %     saveas(gcf, fullfile('fig', ['AutoCovMat_K' num2str(i) '.png']))
% % end
% % 
% % %% reorganise data
% % 
% % % load key data
% % load(fullfile(dir_base,'data','data_DMN_V1-1000_100_TUDA.mat'));
% % load(fullfile(dir_base,'data','T_DMN-1000_100_TUDA.mat'));
% % % load covariates
% % load(fullfile(dir_base,'data','Lat_before_DMN-1000_100_TUDA.mat'));
% % load(fullfile(dir_base,'data','Lat_after_DMN-1000_100_TUDA.mat'));
% % load(fullfile(dir_base,'data','Rsp_moments_DMN-1000_100_TUDA.mat'));
% % load(fullfile(dir_base,'data','CDT_code_DMN-1000_100_TUDA.mat'));
% % 
% % % rewrite data
% % Ts = T; % rewrite T
% % X  = [];
% % Y  = [];
% % T  = [];
% % 
% % % prepare to rewrite data -->  cell to mat
% % 
% % disp('Loading data...')
% % 
% % Sn = length(Ts); Ts_new = cell(20,1); 
% % 
% % for ss = 1:Sn
% %     ttrial = Ts{ss}(1); % assuming the ttrial is the same for all T
% %     xInd = []; n = 0; % counter for T{ss}.
% %     
% %     for i_ntrial = 1:length(Ts{ss})
% %         
% %         % get rid of NaN value; this problem only exits in TUDA analyses
% %         % determine if in this trial, Lat is NaN, if it is, skip
% %         % this trial
% %         if i_ntrial == 1
% %             slice = 1:ttrial;
% %         else
% %             slice = sum(Ts{ss}(1:i_ntrial-1))+1 : sum(Ts{ss}(1:i_ntrial-1)) + ttrial;
% %         end
% %         %--------------------------------------------------------------------------------------
% %         % if isempty(find(isnan(Lat_before{ss}(slice)))) && isempty(find(isnan(Lat_after{ss}(slice))))
% %         % rewrite data to have same length between trials
% %         xInd = [xInd, slice];
% %         n = n+1; % update counter
% %         Ts_new{ss}(n,1) = ttrial;
% %         % end
% %     end
% %     % check inside variables
% %     disp(' ')
% %     disp(['Subj:   ' num2str(ss)])
% %     disp(['Ntrial: ' num2str(length(Ts{ss}))])
% %     disp(['N:      ' num2str(n)])
% %     
% %     data{ss} = data{ss}(xInd,:);
% %     X = [X; data{ss}];
% %     Y = [Y; Lat_before{ss}(xInd,:) Lat_after{ss}(xInd,:) Rsp_moments{ss}(xInd,:) CDT_code{ss}(xInd,:)];
% %     T = [T; Ts_new{ss}];
% %     
% % end % for ss
% % 
% % disp('Data are loaded!')
% % 
% % %% variance check
% % [e_group,e_subj] = explainedvar_PCA(data,Ts_new,options);
% % figure; plot(e_group,'LineWidth',3);
% % save(outputfile, 'e_group', 'e_subj')
% % 
% % %% spectra information
% % % Compute the spectra, at the group level and per subject
% % 
% % N  = 20;
% % Hz = 250;
% % options_mt = struct('Fs',Hz); % Sampling rate - for the 20subj it is 250
% % options_mt.fpass = [1 35];  % band of frequency you're interested in
% % options_mt.tapers = [4 7]; % taper specification - leave it with default values
% % options_mt.p = 0; %0.01; % interval of confidence
% % options_mt.win = 2 * Hz; % multitaper window
% % options_mt.to_do = [1 0]; % turn off pdc
% % options_mt.order = 0;
% % options_mt.embeddedlags = -7:7;
% % 
% % % average
% % fitmt = hmmspectramt(X,T, Gamma,options_mt);
% % save(outputfile, 'fitmt', '-append')
% % %save(outputfile, 'fitmt')
% % 
% % % per subject
% % fitmt_subj = cell(N,1);
% % d = length(options_mt.embeddedlags) - 1;
% % acc = 0;
% % 
% % for n=1:N
% %     gamma = Gamma(acc + (1:(sum(Ts_new{n})-length(Ts_new{n})*d)),:);
% %     acc = acc + size(gamma,1);
% %     fitmt_subj{n} = hmmspectramt(data{n},Ts_new{n},gamma,options_mt);
% %     
% %     fitmt_subj{n}.state = rmfield(fitmt_subj{n}.state,'ipsd');
% %     fitmt_subj{n}.state = rmfield(fitmt_subj{n}.state,'pcoh');
% %     fitmt_subj{n}.state = rmfield(fitmt_subj{n}.state,'phase');
% %     disp(['Subject ' num2str(n)])
% % end
% % save(outputfile,'fitmt_subj','fitmt','-append')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%================================================================
% %% Do an automatic spectral factorisation to find spectrally-defined networks
% % (see paper for more info)
% % Note: manual inspection of the spectral modes (i.e. what is contained in sp_profiles_wb
% % and sp_profiles_4b) is **strongly** recommended. This is an algorithmic
% % solution and there is no theoretical guarantee of getting a sensible
% % result.
% 
% % Get the three bands depicted in the paper (the 4th is essentially capturing noise)
% options_fact = struct();
% options_fact.Ncomp = 4;
% options_fact.Base = 'coh';
% [fitmt_group_fact_4b,sp_profiles_4b,fitmt_subj_fact_4b] = spectdecompose(fitmt_subj,options_fact);
% save(outputfile,'fitmt_subj_fact_4b','fitmt_group_fact_4b','sp_profiles_4b','-append')
% 
% % Get the Four bands depicted in the paper (the 5th is essentially capturing noise)
% options_fact = struct();
% options_fact.Ncomp = 5;
% options_fact.Base = 'coh';
% [fitmt_group_fact_5b,sp_profiles_5b,fitmt_subj_fact_5b] = spectdecompose(fitmt_subj,options_fact);
% save(outputfile,'fitmt_subj_fact_5b','fitmt_group_fact_5b','sp_profiles_5b','-append')
% 
% % Get the Four bands depicted in the paper (the 6th is essentially capturing noise)
% options_fact = struct();
% options_fact.Ncomp = 6;
% options_fact.Base = 'coh';
% [fitmt_group_fact_6b,sp_profiles_6b,fitmt_subj_fact_6b] = spectdecompose(fitmt_subj,options_fact);
% save(outputfile,'fitmt_subj_fact_6b','fitmt_group_fact_6b','sp_profiles_6b','-append')
% 
% % Get the wideband maps (the second is capturing noise)
% options_fact.Ncomp = 2;
% [fitmt_group_fact_wb,sp_profiles_wb,fitmt_subj_fact_wb] = spectdecompose(fitmt_subj,options_fact);
% save(outputfile,'fitmt_subj_fact_wb','fitmt_group_fact_wb','sp_profiles_wb','-append')
% 
% % check if the spectral profiles make sense, if not you might like to repeat
% figure;
% subplot(1,3,1); plot(sp_profiles_4b,'LineWidth',2)
% xlim([0 70])
% subplot(1,3,2); plot(sp_profiles_5b,'LineWidth',2)
% xlim([0 70])
% subplot(1,3,3); plot(sp_profiles_wb,'LineWidth',2)
% xlim([0 70])

% ------ figure -------
load /media/dian/D/data/Binocular_Rivalry/HMM_analyses/output_test/embedded_lag/DMN_V1_1000_100_K4_pca16/run1_analyses/hmm_analysis.mat
addpath(genpath('/media/dian/D/data/scripts/Github/hex_and_rgb_v1.1.1'))
colors = hex2rgb(['798E87', 'C27D38', 'CCC591', '29211F']);

figure;
hax = axes;
xlim([0 40])
hold on

for fqb = 1:4
    plot(sp_profiles_5b(:,fqb), 'Color',colors(fqb,:));
    col = sp_profiles_5b(:,fqb);
    x = find(col == max(col));
    %----- find full width at the half maximum
    hc = max(col)/2;
    diff_hc = abs(col-hc);
    [~, fwhx_ind] = sort(diff_hc);
  %  fwhx = fwhx_ind(1:4);
    fwhm = []; n1=0; n2=0;
    
    for ifx = 1:length(fwhx_ind)
        fx = fwhx_ind(ifx);
        if fx < x && n1==0
            fwhm(1) = fx;
            n1 = 1+n1;
        elseif fx > x && n2==0
            fwhm(2) = fx;
            n2 = n2+1;
        end
    end     
    
    disp(['No. ' num2str(fqb) ' frequency:'])
    disp(['Max frequency = ' num2str(x)])
    disp(['FWHM = [' num2str(fwhm(1)) ' ' num2str(fwhm(2)) ']'])
    disp(' ')
    
    h = area(sp_profiles_5b(:,fqb));
    h.FaceColor = colors(fqb,:) ;
    h.FaceAlpha = 0.4 ;
    line([x x],get(hax,'YLim'),'Color',colors(fqb,:), 'LineWidth',2,'LineStyle','--')
end
hold off
xlabel('Frequency')
ylabel('Power Spectral Density')



%% Do statistical testing on the spectral information
K = options.K; G = 5; N = 20;
within_state = 0;
relative = 0;

fitmt_subj_fact_1d = cell(N,1);
kstats_5b = cell(G-1,1);

for g = 1:G
    for n = 1:N
        fitmt_subj_fact_1d{n} = struct();
        fitmt_subj_fact_1d{n}.state = struct();
        
        for k = 1:K % we don't care about the second component
            fitmt_subj_fact_1d{n}.state(k).psd = fitmt_subj_fact_5b{n}.state(k).psd(g,:,:);
            fitmt_subj_fact_1d{n}.state(k).coh = fitmt_subj_fact_5b{n}.state(k).coh(g,:,:);
        end
        
    end
    
    disp(['Processing for the No. ' num2str(g) ' frequency band...'])
    
    tests_spectra = specttest(fitmt_subj_fact_1d,5000,within_state,relative);
    significant_spectra = spectsignificance(tests_spectra,0.01);
    
    kstats_5b{g}.tests_spectra = tests_spectra;
    kstats_5b{g}.significant_spectra = significant_spectra;
end

save(outputfile,'kstats_5b', '-append')


%------ matrix figures ------
%i_fq = 1;
%k=1;
frq_bands = {'Theta', 'Alpha-1', 'Beta-1', 'Beta-2'};
close all
for i_fq = 1: G-1
    for k=1:K
        
        %psd = squeeze(kstats_5b{freqband}.tests_spectra.higher_corr.state(k).psd);
        %coh = squeeze(kstats_5b{freqband}.tests_spectra.higher_corr.state(k).coh);
        
        psd = squeeze(fitmt_group_fact_5b.state(k).psd(i_fq,:,:));
        coh = squeeze(fitmt_group_fact_5b.state(k).coh(i_fq,:,:));
        
        matx = ones(size(psd));
        matx(psd~=0) = psd(psd~=0);
        matx(coh~=1) = coh(coh~=1);
        
        % ---- plot ----------
        figure; imagesc(matx); hcb=colorbar; colormap(pink);
        %title(hcb, 'P value')
        set(gca,'XTickLabel',ROIs)
        set(gca,'YTickLabel',ROIs)
        tlt = ['K = ' num2str(k), ', Frequency band = ' frq_bands{i_fq}];
        title(tlt)
        
        % 1. find significant position for the upper tail
        psd_sig = squeeze(kstats_5b{i_fq}.significant_spectra.higher_corr.state(k).psd);
        coh_sig = squeeze(kstats_5b{i_fq}.significant_spectra.higher_corr.state(k).coh);
        [r1,c1] = find(psd_sig==1);
        [r2,c2] = find(coh_sig==1);
        rc = [r1 c1; r2 c2];
        
        for ir = 1:size(rc,1)
            pos = rc(ir,:)-0.5;
            rectangle('Position', [pos,1,1], 'EdgeColor','r', 'LineWidth',2)
        end
        
        % 2. find significant position for the lower tail
        psd_sig = squeeze(kstats_5b{i_fq}.significant_spectra.lower_corr.state(k).psd);
        coh_sig = squeeze(kstats_5b{i_fq}.significant_spectra.lower_corr.state(k).coh);
        [r1,c1] = find(psd_sig==1);
        [r2,c2] = find(coh_sig==1);
        rc = [r1 c1; r2 c2];
        
        for ir = 1:size(rc,1)
            pos = rc(ir,:)-0.5;
            rectangle('Position', [pos,1,1], 'EdgeColor','b', 'LineWidth',2)
        end
        
        img_name = ['spectra_k' num2str(k) '_' frq_bands{i_fq} '_xStates.png'];
        saveas(gcf, fullfile('fig',img_name))
    end
end



