%======================== BR GLM0 Extract Covriates =====================
clear % clean up current workspace
%% Define processing steps
glm0       = 1;
estimation = 1;
contrast   = 1;
VOI_covar  = 1;

%% Root Directory
dir_base    = '/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/subs'; %subs/S01/MRI/fMRI_preproc/scan1';
rp_dir_base = '/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/rivalry_fMRI_head_motion_files';%/S01/MRI/fMRI_preproc_rpfiles/scan1'
func        = 'MRI/fMRI_preproc'; % functional scan base
struc       = 'MRI/structural';  % structural scan base

loopin = 1:20;

file = {};

%% Initiating SPM
%spm fmri
%spm('defaults','fmri');
spm_jobman('initcfg');

%% Define subject names
sublist = dir(fullfile(dir_base,'S*'));

%--------------loop base---------------------
% Model Specification: GLM0
if glm0==1    % Model Specification: GLM0
    clear matlabbatch
    for ss = loopin
        sub      = sublist(ss).name;
        subFuncDir  = fullfile(dir_base, sub, func);
        
        scanlist = dir(fullfile(subFuncDir,'scan*'));
        
        for sc = 1:length(scanlist)
            
            scan = scanlist(sc).name;
            FuncDir = fullfile(subFuncDir, scan);
            GLM0_dir = fullfile(FuncDir, 'model', 'GLM0');
            if exist(fullfile(GLM0_dir,'SPM.mat'),'file')==0
                rp_dir = dir(fullfile(rp_dir_base, sub, 'MRI/fMRI_preproc_rpfiles', scan, 'rp*.txt'));
                rp_file = fullfile(rp_dir_base, sub, 'MRI/fMRI_preproc_rpfiles', scan, rp_dir.name);
                rp_outlier_mask = fullfile(rp_dir_base, sub, 'MRI/fMRI_preproc_rpfiles', scan, [sub '_' scan '_outlier_mask.mat']);
                
                rp = load(rp_file);
                rp_outlier = load(rp_outlier_mask);
                rp_outlier.mask(1:5) = 0; % censor the first 5 dummy scans
                
                mkdir(GLM0_dir)
                % select functional scans once for all
                f     = spm_select('List', FuncDir, '^swraf.*\.nii$');
                files = cellstr([repmat([FuncDir '/'],size(f,1),1) f repmat(',1',size(f,1),1)]);
                %movement = load(fullfile(dir_base,sub,sub,'MNINonLinear/Results/tfMRI_WM_LR/Movement_Regressors.txt'));
                % Specify the output directory where the SPM.mat file goes.
                matlabbatch{1}.spm.stats.fmri_spec.dir = {GLM0_dir};
                % Scanner info.
                matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
                matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2.2; % remember to specify
                matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
                matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
                % Pick up the scans. There are 12 sessions for each person
                matlabbatch{1}.spm.stats.fmri_spec.sess.scans = files;
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
                % regressor 1
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name = 'images'; % Have to specify here, to regress with the main effect of the whole images.
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val = 1:length(files); % How many scans
                % regressor 2
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name = 'movement_outlier';
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val = rp_outlier.mask; % outlier mask
                %
                matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {rp_file}; % six movement regressors
                matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
                % Specify experiment parameters
                matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
                matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
                matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
                matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
                matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
                matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
                matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
                
                spm_jobman('run',matlabbatch);
                disp('GLM0 specification completed.');
            end
        end
    end
end
%=================== Estimation, CSF/WM extraction ==================
if estimation ==1
    
    for ss = loopin
        clear matlabbatch
        sub      = sublist(ss).name;
        subFuncDir  = fullfile(dir_base, sub, func);
        
        scanlist = dir(fullfile(subFuncDir,'scan*'));
        
        for sc = 1:length(scanlist)
            
            scan = scanlist(sc).name;
            FuncDir = fullfile(subFuncDir, scan);
            GLM0_dir = fullfile(FuncDir, 'model', 'GLM0');
            if exist(fullfile(GLM0_dir,'mask.nii'),'file')==0
                % Model estimation
                disp(['Estimate GLM0 for Subject: ', sub '-' scan]);
                matlabbatch{1}.spm.stats.fmri_est.spmmat = {[GLM0_dir,'/SPM.mat']}; %Select the SPM.mat
                matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
                matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
                % save batch file for review
                % run batch
                disp('Begin to run...');
                
                spm_jobman('run',matlabbatch);
                disp('GLM0 estimation job completed.');
            end
        end
    end
end

%% Contrast
if contrast == 1
    for ss = loopin
        
        sub      = sublist(ss).name;
        subFuncDir  = fullfile(dir_base, sub, func);
        
        scanlist = dir(fullfile(subFuncDir,'scan*'));
        
        for sc = 1:length(scanlist)
            
            scan = scanlist(sc).name;
            FuncDir = fullfile(subFuncDir, scan);
            GLM0_dir = fullfile(FuncDir, 'model', 'GLM0');
            
            if exist(fullfile(GLM0_dir,'ess_0001.nii'),'file')==0
                clear matlabbatch
                disp(['RUNNING Contrast for subject ' sub '-' scan]);
                % define directories in the session
                matlabbatch{1}.spm.stats.con.spmmat = {[GLM0_dir,'/SPM.mat']};
                matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'all';
                matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights = eye(8);
                matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
                matlabbatch{1}.spm.stats.con.delete = 1;
                % run batch
                disp('Begin to run...');
                spm_jobman('run',matlabbatch);
                disp('''F'' contrast job completed.');
                clear matlabbatch
            end
        end
        
    end
end


%% Extract timeseries of VOIs
if VOI_covar == 1
    for ss = loopin
        clear matlabbatch
        sub      = sublist(ss).name;
        subFuncDir  = fullfile(dir_base, sub, func);
        
        scanlist = dir(fullfile(subFuncDir,'scan*'));
        
        for sc = 1:length(scanlist)
            
            scan = scanlist(sc).name;
            FuncDir = fullfile(subFuncDir, scan);
            GLM0_dir = fullfile(FuncDir, 'model', 'GLM0');
            
            clear matlabbatch
            s{1,1} = '/applications/spm/spm12_7219/tpm/TPM.nii,2';
            s{2,1} = '/applications/spm/spm12_7219/tpm/TPM.nii,3';
            VOIs = {'WM', 'CSF'};
            
            for i = 1:length(VOIs)
                voi = VOIs{i};
                
                if exist(fullfile(GLM0_dir,['VOI_' voi '_1.mat']),'file')==0
                    
                    clear matlabbatch Y
                    disp(['Extract ' VOIs{i} ' Timeseries  for subject ', sub '-' scan '.']);
                    matlabbatch{1}.spm.util.voi.spmmat = {[GLM0_dir '/SPM.mat']};
                    matlabbatch{1}.spm.util.voi.adjust = NaN;
                    matlabbatch{1}.spm.util.voi.session = 1;
                    matlabbatch{1}.spm.util.voi.name = voi;
                    matlabbatch{1}.spm.util.voi.roi{1}.mask.image = s(i,1);
                    matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.7; % avoid mixed tissues because of the smoothing (CésarCaballero-GaudesaRichard C.Reynoldsb., 2016)
                    matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {fullfile(GLM0_dir, 'mask.nii,1')};
                    matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
                    matlabbatch{1}.spm.util.voi.expression = 'i1&i2';
                    % run batch
                    spm_jobman('run',matlabbatch);
                    disp([VOIs{i} ' extraction job completed.']);
                    clear matlabbatch
                end
            end
            clear i
            if exist(fullfile(FuncDir,'allconfounds.txt'),'file')==0
            cd(FuncDir)
            disp('Writing confounds..');
            allconfounds = [];
            load(fullfile(GLM0_dir, 'VOI_WM_1.mat'))
            allconfounds(:,1) = Y;
            clear Y
            load(fullfile(GLM0_dir, 'VOI_CSF_1.mat'))
            allconfounds(:,2) = Y;
            
            % load head movement
            rp_dir = dir(fullfile(rp_dir_base, sub, 'MRI/fMRI_preproc_rpfiles', scan, 'rp*.txt'));
            rp_file = fullfile(rp_dir_base, sub, 'MRI/fMRI_preproc_rpfiles', scan, rp_dir.name);
            rp = load(rp_file);
            
            allconfounds = [allconfounds rp];
            % save file
            dlmwrite(fullfile(GLM0_dir, 'allconfounds.txt'),allconfounds,'delimiter',' ')
            disp('Done.')
            end
        end
    end
end




