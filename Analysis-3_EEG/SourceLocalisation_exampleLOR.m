%-----------------------------------------------------------------------
% BR source localisation using batch script
% spm SPM - SPM12 (7219); 11 Aug 2019, by Dian Lu
% 1. headmodel specification (MRI coord manually selected, QCed and typed) + forward model
% 2. source inverse
%-----------------------------------------------------------------------
tic
clear
rmpath(genpath('/applications/spm/spm12_6906')) % VERY important! Old version fails.
addpath('/home/dl577/spm12')
addpath('/data/dl577/scripts/matlab/self-defined_functions/')
spm('defaults', 'EEG');
% load info of condition-session correspondence for later use
load('/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/info/replay_inst_ind.mat'); % a cell
load('/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/info/replay_smooth_ind.mat');
load('/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/info/rival_ind.mat');
%% what analyses?
headmodel = 1;
inverse   = 1;

%% predefine MRI coord
cord(:,:,1) = [-75.6 23.0 -8.8; 73.7 6.8 -21.9; -3.4 129.3 53.7]; % subj1=[TP9; TP10; Fpz]
cord(:,:,2) = [-73.4 8.4 -19.6; 80.3 11.5 -22.8; 3.9 151.3 30.2];
cord(:,:,3) = [-67.0 7.0 -39.2; 79.7 17.6 -35.0; -5.0 144.8 43.1];
cord(:,:,4) = [-75.0 -1.0 -35.0; 73.0 4.0 -29.0; 2.5 152.5 24.5];
cord(:,:,5) = [-79.3 23.3 -37.0; 73.1 23.3 -23.2; -5.2 157.1 48.7];
cord(:,:,6) = [-76.3 14.9 -17.7; 74.9 7.5 -20.8; 4.5 147.9 38.0];
cord(:,:,7) = [-70.3 23.0 -44.3; 84.7 17.7 -51.6; 8.8 133.5 36.9];
cord(:,:,8) = [-68.6 -4.5 -37.8; 80.9 11.8 -37.9; -4.5 119.6 44.5];
cord(:,:,9) = [-69.6 10.0 -40.0; 67.9 20.5 -42.1; -1.4 155.4 23.9];
cord(:,:,10) = [-79.7 -53.0 -33.1; 72.2 -38.1 -45.9; 0.0 67.9 48.6];
cord(:,:,11) = [-77.2 16.4 -33.9; 80.8 26.5 -28.9; -3.3 133.5 64.2];
cord(:,:,12) = [-77.2 14.0 -39.3; 70.7 9.8 -26.7; 0.4 150.3 27.8];
cord(:,:,13) = [-63.6 -3.7 -59.9; 78.1 -6.1 -49.7; 1.1 127.3 13.1];
cord(:,:,14) = [-78.7 -46.2 -33.9; 77.9 -38.9 -37.0; -3.0 73.1 45.0];
cord(:,:,15) = [-70.2 10.1 -20.9; 65.1 1.1 -8.4; 0.6 139.6 38.4];
cord(:,:,16) = [-71.3 -40.8 -20.9; 67.8 -44.9 -25.9; 0.7 85.3 21.1];
cord(:,:,17) = [-68.8 23.5 -33.0; 72.4 22.4 -39.6; 5.1 154.2 51.7];
cord(:,:,18) = [-72.3 7.5 -23.2; 79.4 0.3 -31.4; 21.0 118.9 40.2];
cord(:,:,19) = [-76.3 -28.0 -40.3; 64.8 -46.2 -33.9; 6.0 81.0 29.2];
cord(:,:,20) = [-77.9 -51.4 -15.4; 76.1 -43.2 -21.6; 0.6 79.9 26.0];

%% define base directory
dir_base = '/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/EEG/subs-eeg';
GLM2_base = '/data/dl577/BR/activation_AR';
files = {'/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/EEG/SPM_objects/CS01_scan1.mat'};

%% loop for each epoch files in sub dirs.
for i = 1%472:571%:length(files)
    file = files{i};
    fn = strsplit(file,'/');
    sub = 'S01';%fn{end-1};
    %  find Sub index
    noS = sub(2:3); noS(find(noS=='0')) = [];
    inS = str2num(noS);
    matfile = fn{end};
    subdir = fullfile(dir_base, sub);
    mri_file = fullfile(subdir, 'MRI/structural', [sub, '_MRI.nii,1']);
    
    % prior in-use depends on the condition of the current session
    strs = strsplit(matfile,'scan');
    epoch_prefix = strs{1}(1:2);
    scan = strsplit(strs{2}, '.mat'); inScan = str2num(scan{1});
    
    RplSM = replay_smooth_ind{inS};
    RplIT = replay_inst_ind{inS};
    BR    = rival_ind{inS};
    
    if ismember(inScan, BR) && strcmp(epoch_prefix,'ed')
        model = 'SSBR-baseline-grnrd';
        clear pfile
        pfile = dir(fullfile(GLM2_base, model, 'prior*.nii'));
        pfile = fullfile(GLM2_base, model, pfile.name);
        
    elseif ismember(inScan, BR) && strcmp(epoch_prefix,'eg')
        model = 'BR-baseline-green';
        clear pfile
        pfile = dir(fullfile(GLM2_base, model, 'prior*.nii'));
        pfile = fullfile(GLM2_base, model, pfile.name);
        
    elseif ismember(inScan, BR) && strcmp(epoch_prefix,'er')
        model = 'BR-baseline-red';
        clear pfile
        pfile = dir(fullfile(GLM2_base, model, 'prior*.nii'));
        pfile = fullfile(GLM2_base, model, pfile.name);
        
    elseif ismember(inScan, BR) && strcmp(epoch_prefix,'eu')
        model = 'BR-baseline-mix';
        clear pfile
        pfile = dir(fullfile(GLM2_base, model, 'prior*.nii'));
        pfile = fullfile(GLM2_base, model, pfile.name);
        %-----------------------------------------------------
    elseif ismember(inScan, RplSM) && strcmp(epoch_prefix,'ed')
        model = 'SSRpl-baseline-grnrd';
        pfile = dir(fullfile(GLM2_base, model, 'prior*.nii'));
        pfile = fullfile(GLM2_base, model, pfile.name);
        
    elseif ismember(inScan, RplSM) && strcmp(epoch_prefix,'eg')
        model = 'RplSM-baseline-green';
        clear pfile
        pfile = dir(fullfile(GLM2_base, model, 'prior*.nii'));
        pfile = fullfile(GLM2_base, model, pfile.name);
        
    elseif ismember(inScan, RplSM) && strcmp(epoch_prefix,'er')
        model = 'RplSM-baseline-red';
        clear pfile
        pfile = dir(fullfile(GLM2_base, model, 'prior*.nii'));
        pfile = fullfile(GLM2_base, model, pfile.name);
        
    elseif ismember(inScan, RplSM) && strcmp(epoch_prefix,'eu')
        model = 'RplSM-baseline-mix';
        clear pfile
        pfile = dir(fullfile(GLM2_base, model, 'prior*.nii'));
        pfile = fullfile(GLM2_base, model, pfile.name);
        %-----------------------------------------------------
    elseif ismember(inScan, RplIT)
        disp(['Skipping ' matfile '...'])
        continue % at this point, skip Rpl Instantaneous condition
    end
    
   % cd(subdir) % subj folder
    
    %% Head model specification + forward model
    if headmodel
        clear matlabbatch
        matlabbatch{1}.spm.meeg.source.headmodel.D = {file};
        matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
        matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
        matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.mri = {mri_file};
        matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'TP9';
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = cord(1,:,inS);
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'TP10';
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = cord(2,:,inS);
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'Fpz';
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = cord(3,:,inS);
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 1;
        matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = '3-Shell Sphere';
        matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
        % run batch
        spm_jobman('run', matlabbatch);
        disp(['Head-model specified and foward model solved for ' matfile '.'])
    end
    
    %% source inversion
    if inverse
        %diary(['log_source_localisation_' sub '.txt'])
       % diary on
        % write matlatbbatch
        clear matlabbatch
        matlabbatch{1}.spm.meeg.source.invert.D = {file};
        matlabbatch{1}.spm.meeg.source.invert.val = 1;
        matlabbatch{1}.spm.meeg.source.invert.whatconditions.all = 1;
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.invtype = 'LOR';
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.woi = [-Inf Inf];
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.foi = [1 35]; % adjust for different study
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.hanning = 1;
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.priors.priorsmask = {pfile};
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.priors.space = 1; % MNI space
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.locs = zeros(0, 3);
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.radius = 32;
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.mask = {''};
        matlabbatch{1}.spm.meeg.source.invert.modality = {'All'};
        % run batch
        spm_jobman('run', matlabbatch);
        disp(['Source inversion done for ' matfile '.'])
        diary off
    end
end
toc
t = toc