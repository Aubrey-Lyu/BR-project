%% Art_global repair method to deal with head movement
% rmpath(genpath('/applications/spm/spm12_6906'))
% addpath('/home/dl577/spm12')
% addpath('/home/dl577/spm12/toolbox/ArtRepair')
% function art_global(Images,RealignmentFile,HeadMaskType, RepairType)
% HeadMaskType  = 1 for SPM mask, = 4 for Automask
%    RepairType = 1 for ArtifactRepair alone (0.5 movement and add margin).
%               = 2 for Movement Adjusted images  (0.5 movement, no margin)
%               = 0 No repairs are done, bad scans are found.
%                   Listed in art_suspects.txt for motion adjustment.
%% Root Directory
dir_base    = '/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/subs'; %subs/S01/MRI/fMRI_preproc/scan1';
rp_dir_base = '/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/rivalry_fMRI_head_motion_files';%/S01/MRI/fMRI_preproc_rpfiles/scan1'
func        = 'MRI/fMRI_preproc'; % functional scan base
struc       = 'MRI/structural';  % structural scan base

%% Define subject names
sublist = dir(fullfile(dir_base,'S*'));
loopin = 1:20;
%--------------loop base---------------------
for ss = loopin
    sub      = sublist(ss).name;
    subFuncDir  = fullfile(dir_base, sub, func);
    
    scanlist = dir(fullfile(subFuncDir,'scan*'));
    
    for sc = 1:length(scanlist)
        
        scan = ['scan' num2str(sc)];
        FuncDir = fullfile(subFuncDir, scan);
        
        clear rp f
        %if exist(fullfile(FuncDir,'ArtifactMask.nii'),'file')==0
            rp_dir = dir(fullfile(rp_dir_base, sub, 'MRI/fMRI_preproc_rpfiles', scan, 'rp*.txt'));
            rp_file = fullfile(rp_dir_base, sub, 'MRI/fMRI_preproc_rpfiles', scan, rp_dir.name);
            rp_outlier_mask = fullfile(rp_dir_base, sub, 'MRI/fMRI_preproc_rpfiles', scan, [sub '_' scan '_outlier_mask.mat']);
            % select functional scans once for all
            file = {};
            f     = spm_select('List', FuncDir, '^swraf.*\.nii$');
            %files = cellstr([repmat([FuncDir '/'],size(f,1),1) f repmat(',1',size(f,1),1)]);
            files = [repmat([FuncDir '/'],size(f,1),1) f repmat(',1',size(f,1),1)];
            disp([sub '-' scan ' has ' num2str(size(files,1)) ' volumes...'])
            % apply art_global
           % art_global(files, rp_file, 4, 2)
           % disp(['Art repaired fast-motion file for ' sub '-' scan '.'])
        %end
    end
end



