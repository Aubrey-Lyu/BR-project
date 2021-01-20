%----------Second level analysis: BR activation ---------------------------
% Job saved on 06-Aug-2019 12:04:01 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (6906)
% cfg_basicio BasicIO - Unknown
% some of the scans are art-repaired; and extra contrasts were added.
%% Scheme:
% 1. 1-10 one-sample T for main effect of baselines and 1stlevel-contrasts
% 2. 11-13 paired T for contrasting between 1st-level baselines
% 3. 14 anova for testing the BR(r&g)-Replay(r&g) difference
% 4. 15 one-sample T test for main effect of stable percepts (r&g)
%-----------------------------------------------------------------------
clear
addpath ~/scripts/spm_scripts/BR
%% statistics to do
onesampleT  = 0;
pairT       = 0;
anova       = 0;
onesampleT2 = 1; % for the combination of green and red percepts
estimate    = 0;

%% load contrast names
load('/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/subs/S01/MRI/fMRI_preproc/model/GLM1_bhv/SPM.mat')
xCon = {SPM.xCon.name};

%% load info of subs and scans that need to be repaired
aT = import_improved_scans('/data/dl577/BR/info/improved_scans.txt');
aT = aT(strcmp(aT.Region, 'total'),:); % select only "total" region comparison
%% define variables/folders
dir_base    = '/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/subs'; %subs/S01/MRI/fMRI_preproc/scan1';
output_base = '/data/dl577/BR/activation_AR';
mkdir(output_base);
sublist     = dir(fullfile(dir_base,'S*'));
func        = 'MRI/fMRI_preproc'; % functional scan base

loopin = 1:20;

%% Initiating SPM
%spm fmri
%spm('defaults','fmri');
spm_jobman('initcfg');

%% One-sample T test
if onesampleT == 1
    for c = 1:length(xCon)
        cname = xCon{c};
        output_folder = fullfile(output_base, cname);
        mkdir(output_folder);
        
        clear matlabbatch
        matlabbatch{1}.spm.stats.factorial_design.dir = {output_folder};
        % prepare scans
        scans = {};
        for ss = loopin
            sub      = sublist(ss).name;
            if isempty(find(aT.Sub==ss & aT.image==c)) % if not in the repair list
                GLM1_dir = fullfile(dir_base, sub, func, 'model', 'GLM1_bhv');
            else
                GLM1_dir = fullfile(dir_base, sub, func, 'model', 'GLM1_bhv_ArtRepair');
            end
            scans{end+1,1} = fullfile(GLM1_dir, ['con_' sprintf('%04d',c) '.nii,1']);
        end
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = scans;
        % other parameters
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        % run
        spm_jobman('run',matlabbatch);
        disp(['One sample T done for the no. ' num2str(c) ' contrast: ' cname '.'])
    end
end

%% Paired T test
if pairT == 1
    pCon = xCon(1:3);
    pCon = strrep(pCon,'BR-baseline','SSBR-RplSM'); % define con names
    for c = 1:length(pCon)
        cname = pCon{c};
        output_folder = fullfile(output_base, cname);
        mkdir(output_folder);
        
        clear matlabbatch
        matlabbatch{1}.spm.stats.factorial_design.dir = {output_folder};
        
        % prepare scans
        scans = {};
        for ss = loopin
            sub      = sublist(ss).name;
            % for the 1st con image
            if isempty(find(aT.Sub==ss & aT.image==c)) % if not in the repair list
                GLM1_dir1 = fullfile(dir_base, sub, func, 'model', 'GLM1_bhv');
            else
                GLM1_dir1 = fullfile(dir_base, sub, func, 'model', 'GLM1_bhv_ArtRepair');
            end
            % for the 2nd con image
            if isempty(find(aT.Sub==ss & aT.image==c+3)) % if not in the repair list
                GLM1_dir2 = fullfile(dir_base, sub, func, 'model', 'GLM1_bhv');
            else
                GLM1_dir2 = fullfile(dir_base, sub, func, 'model', 'GLM1_bhv_ArtRepair');
            end
            % select scans
            scans = {fullfile(GLM1_dir1, ['con_' sprintf('%04d',c) '.nii,1']);... % 1,2,3
                fullfile(GLM1_dir2, ['con_' sprintf('%04d',c+3) '.nii,1'])}; % 4,5,6
            matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(ss).scans = scans; % ss (no. of sbjs) pairs
        end
        
        matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0;
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        % run
        spm_jobman('run',matlabbatch);
        disp(['Paired T done for the contrast: ' cname '.'])
    end
end

%% anova
if anova == 1
    cname = 'SSBR-RplSM-grnrd';
    output_folder = fullfile(output_base, cname);
    mkdir(output_folder);
    %
    clear matlabbatch
    matlabbatch{1}.spm.stats.factorial_design.dir = {output_folder};
    scans1 = {}; scans2 = {}; scans3 = {}; scans4 = {};
    for ss = loopin
        sub      = sublist(ss).name;
        % for con_0001.nii
        if isempty(find(aT.Sub==ss & aT.image==1)) % if not in the repair list
            GLM1_dir1 = fullfile(dir_base, sub, func, 'model', 'GLM1_bhv');
        else
            GLM1_dir1 = fullfile(dir_base, sub, func, 'model', 'GLM1_bhv_ArtRepair');
        end
        % for con_0002.nii
        if isempty(find(aT.Sub==ss & aT.image==2)) % if not in the repair list
            GLM1_dir2 = fullfile(dir_base, sub, func, 'model', 'GLM1_bhv');
        else
            GLM1_dir2 = fullfile(dir_base, sub, func, 'model', 'GLM1_bhv_ArtRepair');
        end
        % for con_0004.nii
        if isempty(find(aT.Sub==ss & aT.image==4)) % if not in the repair list
            GLM1_dir3 = fullfile(dir_base, sub, func, 'model', 'GLM1_bhv');
        else
            GLM1_dir3 = fullfile(dir_base, sub, func, 'model', 'GLM1_bhv_ArtRepair');
        end
        % for con_0005.nii
        if isempty(find(aT.Sub==ss & aT.image==5)) % if not in the repair list
            GLM1_dir4 = fullfile(dir_base, sub, func, 'model', 'GLM1_bhv');
        else
            GLM1_dir4 = fullfile(dir_base, sub, func, 'model', 'GLM1_bhv_ArtRepair');
        end
        
        scans1{end+1,1} = fullfile(GLM1_dir1, 'con_0001.nii,1');
        scans2{end+1,1} = fullfile(GLM1_dir2, 'con_0002.nii,1');
        scans3{end+1,1} = fullfile(GLM1_dir3, 'con_0004.nii,1');
        scans4{end+1,1} = fullfile(GLM1_dir4, 'con_0005.nii,1');
        
    end
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(1).scans = scans1;
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(2).scans = scans2;
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(3).scans = scans3;
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(4).scans = scans4;
    matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    % run
    spm_jobman('run',matlabbatch);
    disp(['ANOVA done for the contrast: ' cname '.'])
end

%% One-sample T for both r and g stable percepts.
if onesampleT2 == 1
   % cnames = {'SSBR-baseline-grnrd', 'SSRpl-baseline-grnrd';...
       % [1 2], [4 5]};
       
       cnames = {'mixed-dominant_RplSM', ;...
       17};
       
        for cn = 1:size(cnames,2)
        cname = cnames{1,cn};
        c_i1 = cnames{2,cn}(1);
        %c_i2 = cnames{2,cn}(2);
        output_folder = fullfile(output_base, cname);
        mkdir(output_folder);
        
        clear matlabbatch
        matlabbatch{1}.spm.stats.factorial_design.dir = {output_folder};
        % prepare scans
        scans = {};
        for ss = loopin
            sub      = sublist(ss).name;
            % for con_0001.nii or con_0004.nii
            if isempty(find(aT.Sub==ss & aT.image==c_i1)) % if not in the repair list
                GLM1_dir1 = fullfile(dir_base, sub, func, 'model', 'GLM1_bhv');
            else
                GLM1_dir1 = fullfile(dir_base, sub, func, 'model', 'GLM1_bhv_ArtRepair');
            end
%             % for con_0002.nii or con_0005.nii
%             if isempty(find(aT.Sub==ss & aT.image==c_i2)) % if not in the repair list
%                 GLM1_dir2 = fullfile(dir_base, sub, func, 'model', 'GLM1_bhv');
%             else
%                 GLM1_dir2 = fullfile(dir_base, sub, func, 'model', 'GLM1_bhv_ArtRepair');
%             end
            scans{end+1,1} = fullfile(GLM1_dir1, ['con_00' num2str(c_i1) '.nii,1']); % green
%             scans{end+1,1} = fullfile(GLM1_dir2, ['con_000' num2str(c_i2) '.nii,1']); % red
        end
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = scans;
        % other parameters
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        % run
        spm_jobman('run',matlabbatch);
        disp(['One sample T done for the contrast: ' cname '.'])
        end
end

%%
if estimate == 1
    addpath ~/scripts/matlab/self-defined_functions
    spm_files = whereis(output_base,1,'SPM.mat'); % find the spm path.
    for i = 1:length(spm_files)
        spmfile = spm_files{i};
        
        clear matlabbatch
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {spmfile}; %Select the SPM.mat of glm1
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        % run batch
        disp('Begin to run...');
        spm_jobman('run',matlabbatch);
        disp(['Estimation done for the No. ' num2str(i) ' SPM out of a total of ' num2str(length(spm_files)) ' files.']);
    end
end