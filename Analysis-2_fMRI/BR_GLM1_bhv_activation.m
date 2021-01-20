%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bhv choice--> activation: BR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
%rmpath(genpath('/applications/spm/spm12_6906'))
%addpath('/home/dl577/spm12')
%% WHAT ANALYSES?
glm1     = 0;
estimate = 0;
contrast = 1;

%% Root Directory
dir_base    = '/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/subs'; %subs/S01/MRI/fMRI_preproc/scan1';
rp_dir_base = '/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/rivalry_fMRI_head_motion_files';%/S01/MRI/fMRI_preproc_rpfiles/scan1'
func        = 'MRI/fMRI_preproc'; % functional scan base
struc       = 'MRI/structural';  % structural scan base

loopin = 1:20;

file = {};
%% load scan-condition information
load('/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/info/replay_smooth_ind.mat') % index of scans/sessions
load('/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/info/rival_ind.mat')
%% Initiating SPM
%spm fmri
%spm('defaults','fmri');
spm_jobman('initcfg');

%% Define subject names
sublist = dir(fullfile(dir_base,'S*'));

%--------------loop base---------------------
% Model Specification: GLM1
if glm1==1
    
    for ss = loopin
        sub      = sublist(ss).name;
        subFuncDir  = fullfile(dir_base, sub, func);
        
        scanlist = dir(fullfile(subFuncDir,'scan*'));
        GLM1_dir = fullfile(subFuncDir, 'model', 'GLM1_bhv');
        mkdir(GLM1_dir);
        cd(GLM1_dir);
        
        %   if exist(fullfile(GLM1_dir,'SPM.mat'),'file')==0 % better exist in case of program crushing
        disp(['Specify Model for Subject: ', sub]);
        clear matlabbatch
        % cfg_basicio BasicIO - Unknown
        %-----------------------------------------------------------------------
        matlabbatch{1}.spm.stats.fmri_spec.dir = {GLM1_dir};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        if ss == 13
            tr = 2;
        else
            tr = 2.2;
        end
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = tr;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        for sc = 1:length(scanlist)
            
            scan = ['scan' num2str(sc)]; % order matters here
            FuncDir = fullfile(subFuncDir, scan);
            
            rp_dir = dir(fullfile(rp_dir_base, sub, 'MRI/fMRI_preproc_rpfiles', scan, 'rp*.txt'));
            rp_file = fullfile(rp_dir_base, sub, 'MRI/fMRI_preproc_rpfiles', scan, rp_dir.name);
            
            % Get ready the functional scans
            
            clear f
            files = {};
            % select functional scans once for all
            f     = spm_select('List', FuncDir, '^swraf.*\.nii$');
            % discard the first 5 scans
            files = cellstr([repmat([FuncDir '/'],size(f,1),1) f repmat(',1',size(f,1),1)]);
            
            
            %% load bhv time point
            for i =1:5
                load(fullfile(dir_base, sub, 'bhv', [sub '_' scan '_resp' num2str(i) '.mat']))
            end
            load(fullfile(dir_base, sub, 'bhv', [sub '_' scan '_time_trial.mat']))
            load(fullfile(dir_base, sub, 'bhv', [sub '_' scan '_time_tr.mat']))
            
            rp1_choice = resp1(:,1);
            rp1_time = resp1(:,2) + time_trial(1);
            rp1_diff = [diff(resp1(:,2)); NaN];
            
            rp2_choice = resp2(:,1);
            rp2_time = resp2(:,2) + time_trial(2);
            rp2_diff = [diff(resp2(:,2)); NaN];
            
            rp3_choice = resp3(:,1);
            rp3_time = resp3(:,2) + time_trial(3);
            rp3_diff = [diff(resp3(:,2)); NaN];
            
            rp4_choice = resp4(:,1);
            rp4_time = resp4(:,2) + time_trial(4);
            rp4_diff = [diff(resp4(:,2)); NaN];
            
            rp5_choice = resp5(:,1);
            rp5_time = resp5(:,2) + time_trial(5);
            rp5_diff = [diff(resp5(:,2)); NaN];
            
            choice = [rp1_choice; rp2_choice; rp3_choice; rp4_choice; rp5_choice];
            
            time = [rp1_time; rp2_time; rp3_time; rp4_time; rp5_time];
            diff_time = [rp1_diff ; rp2_diff; rp3_diff; rp4_diff; rp5_diff];
            % head movement file
            rp=load(rp_file);
            % CSF, WM
            load(fullfile(FuncDir, 'model', 'GLM0', 'VOI_CSF_1.mat'))
            CSF = Y;
            load(fullfile(FuncDir, 'model', 'GLM0', 'VOI_WM_1.mat'))
            WM = Y;
            % remember to discard the first 5 scans
            if length(time_tr) ~= 129
                display(['Caution! The ' sub '-' scan ' has ' num2str(length(time_tr)) ' recording of volumes!'])
                %                      if length(time_tr) < size(files,1)
                %                          files = files(1:length(time_tr));
                %                      end
            end
            T = table(choice, time, diff_time);
            
            %                 % apply to data related to native scans, delete the first
            %                 % 11 to make correspondence with the EEG recording.
            %                 files = files(12:end);
            %                 rp = rp(12:end,:); CSF = CSF(12:end); WM = WM(12:end);
            
            %%  put values into experimental design
            Conditions = {'green', 'red', 'mixed', 'task', 'rest'; % cond(sc).name
                T.time(T.choice==1), T.time(T.choice==2), T.time(T.choice==3),... % cond(sc).onset: green, red, mixed
                time_trial, time_trial+42; % task onset (42s task block), rest onset
                0, 0, 0, 42, 12 ... % cond(sc).dur; spike event
                };
            
            task_onset_scan = round(time_trial/tr);
            half_rest_scan = round(6/tr); % 12s rest before and after a task block
            %                 task_rest_contrast = -ones(1, size(files,1));
            %                 for i = 1:5
            %                     task_rest_contrast(task_onset_scan(i):task_onset_scan(i)+20) = 1;
            %                 end
            % create covariates of different means for different runs
            RUN1 = zeros(1, size(files,1)); RUN2 = zeros(1, size(files,1)); RUN3 = zeros(1, size(files,1));
            RUN4 = zeros(1, size(files,1)); RUN5 = zeros(1, size(files,1));
            RUN1(task_onset_scan(1)-half_rest_scan: task_onset_scan(2)-half_rest_scan) = 1;
            RUN2(task_onset_scan(2)-half_rest_scan: task_onset_scan(3)-half_rest_scan) = 1;
            RUN3(task_onset_scan(3)-half_rest_scan: task_onset_scan(4)-half_rest_scan) = 1;
            RUN4(task_onset_scan(4)-half_rest_scan: task_onset_scan(5)-half_rest_scan) = 1;
            RUN5(task_onset_scan(5)-half_rest_scan: task_onset_scan(5)+ round((42+6)/tr)) = 1;
            
            Regressions = { 'run1','run2','run3', 'run4', 'run5','WM','CSF',...
                'rp1', 'rp2' ,'rp3','rp4','rp5','rp6';
                RUN1, RUN2, RUN3, RUN4, RUN5, WM, CSF,...
                rp(:,1), rp(:,2), rp(:,3), rp(:,4),rp(:,5), rp(:,6)};
            
            %% GLM specification
            matlabbatch{1}.spm.stats.fmri_spec.sess(sc).scans = files;
            for cdt = 1:size(Conditions, 2)
                matlabbatch{1}.spm.stats.fmri_spec.sess(sc).cond(cdt).name = Conditions{1, cdt};
                matlabbatch{1}.spm.stats.fmri_spec.sess(sc).cond(cdt).onset = Conditions{2, cdt};
                matlabbatch{1}.spm.stats.fmri_spec.sess(sc).cond(cdt).duration = Conditions{3, cdt};
                matlabbatch{1}.spm.stats.fmri_spec.sess(sc).cond(cdt).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sc).cond(cdt).pmod = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(sc).cond(cdt).orth = 1;
            end
            % multicondition--> none
            matlabbatch{1}.spm.stats.fmri_spec.sess(sc).multi = {''};
            % regression
            for rg = 1:size(Regressions, 2)
                matlabbatch{1}.spm.stats.fmri_spec.sess(sc).regress(rg).name = Regressions{1, rg};
                matlabbatch{1}.spm.stats.fmri_spec.sess(sc).regress(rg).val = Regressions{2, rg};
            end
            % multiple regressor: specify movement regressors; rp file
            % is modified...so not using it
            matlabbatch{1}.spm.stats.fmri_spec.sess(sc).multi_reg = {''};
            % high-pass filter
            matlabbatch{1}.spm.stats.fmri_spec.sess(sc).hpf = 128;
        end
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        
        % run batch
        disp('Begin to run...');
        spm_jobman('run',matlabbatch);
    end
end
%end

%% Model estimation
if estimate == 1
    for ss = loopin
        sub      = sublist(ss).name;
        subFuncDir  = fullfile(dir_base, sub, func);
        
        scanlist = dir(fullfile(subFuncDir,'scan*'));
        
        clear matlabbatch
        
        GLM1_dir = fullfile(subFuncDir, 'model', 'GLM1_bhv');
        
        disp(['Specify Model for Subject: ', sub]);
        
        if exist(fullfile(GLM1_dir,'mask.nii'),'file')==0 && exist(fullfile(GLM1_dir,'SPM.mat'),'file')==2
            
            disp(['Estimate glm1 for Subject: ', sub]);
            
            matlabbatch{1}.spm.stats.fmri_est.spmmat = {[GLM1_dir,'/SPM.mat']}; %Select the SPM.mat of glm1
            matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            
            % run batch
            disp('Begin to run...');
            spm_jobman('run',matlabbatch);
        end
    end
end

%% T Contrasts for First-level
if contrast == 1
    for ss = loopin
        sub      = sublist(ss).name;
        subFuncDir  = fullfile(dir_base, sub, func);

        GLM1_dir = fullfile(subFuncDir, 'model', 'GLM1_bhv');
        glmfile = fullfile(GLM1_dir,'SPM.mat');
        %% preparing contrast values
        
        clear SPM.mat; load(glmfile)
        % for mixed-dominant (BR)
        contrast1 = zeros(1,length(SPM.Vbeta));
        contrast1(3+(rival_ind{ss}-1)*18)=2; % BR-mixed
        contrast1(1+(rival_ind{ss}-1)*18)=-1; % BR-green
        contrast1(2+(rival_ind{ss}-1)*18)=-1; % BR-red
        % for mixed-dominant (RplSM)
        contrast2 = zeros(1,length(SPM.Vbeta));
        contrast2(3+(replay_smooth_ind{ss}-1)*18)=2; % RplSM-mixed
        contrast2(1+(replay_smooth_ind{ss}-1)*18)=-1; % RplSM-green
        contrast2(2+(replay_smooth_ind{ss}-1)*18)=-1; % RplSM-red
        % for mixed-dominant (BR - RplSM)
        contrast3 = zeros(1,length(SPM.Vbeta));
        contrast3(1+(rival_ind{ss}-1)*18)=-1;
        contrast3(2+(rival_ind{ss}-1)*18)=-1;
        contrast3(3+(rival_ind{ss}-1)*18)=2;
        contrast3(1+(replay_smooth_ind{ss}-1)*18)=1;
        contrast3(2+(replay_smooth_ind{ss}-1)*18)=1;
        contrast3(3+(replay_smooth_ind{ss}-1)*18)=-2;
       
        % make a meta CONTRAST file
        CONTRASTS = {'mixed-dominant_BR', 'mixed-dominant_RplSM', 'mixed-dominant_BR-RplSM';... % contrast names
             contrast1, contrast2, contrast3 % condition index or contrast values
            };

        disp(['Beta contrasts for Subject: ', sub]);
        
        %% write matlabbatch
        if exist(fullfile(GLM1_dir,'ResMS.nii'),'file')==2
            
            clear matlabbatch
            
            matlabbatch{1}.spm.stats.con.spmmat = {glmfile};
            for ctr = 1:size(CONTRASTS,2)

                    matlabbatch{1}.spm.stats.con.consess{ctr}.tcon.name = CONTRASTS{1,ctr};
                    matlabbatch{1}.spm.stats.con.consess{ctr}.tcon.weights =  CONTRASTS{2,ctr};
                    matlabbatch{1}.spm.stats.con.consess{ctr}.tcon.sessrep = 'none';
            end
            matlabbatch{1}.spm.stats.con.delete = 0; % do not delete
        %    
        else
            disp('Estimation seems not completed.')
        end
        
        % run batch
        disp('Begin to run...');
        spm_jobman('run',matlabbatch);

    end
end



