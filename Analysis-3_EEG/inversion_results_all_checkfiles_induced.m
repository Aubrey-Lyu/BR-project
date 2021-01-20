%-----------------------------------------------------------------------
% save inversion results as a image, then do statisical testing
%-----------------------------------------------------------------------
clear

%% what analyses
inv_nii = 1;

%% set up
rmpath(genpath('/applications/spm/spm12_6906'))
addpath('/home/dl577/spm12')
addpath ~/scripts/matlab/self-defined_functions/
spm('defaults', 'EEG');
spm_jobman('initcfg');

% load('/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/info/rival_ind.mat');
% load('/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/info/replay_smooth_ind.mat'); % index of scans/sessions
dir_base = '/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/EEG/subs-eeg';
percepts = { 's'; 'u'};
tic
% keep diary
%diary /lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/info/log_inversion_write_nii.txt
%diary on

%% save result for source reconstruction
if inv_nii
    Inv = [2 7 9 11 13 15 17]; % inversion model index
    Woi = [-Inf Inf; -0.2*1000 0.2*1000]; % in ms
    Foi = [1 35; 1 3.5; 4 7; 7.5 9.5; 10 12; 13 23; 24 34];
    pars = struct('Inv', Inv, 'Foi', Foi); % these two pars are bound together
    %% loop though these
    loopin_w = 1:size(Woi,1);
    loopin_i = 1:length(pars.Inv);
    loopin_s = 1:20;
    %% set configurations
    for w = loopin_w
        for i = loopin_i
            disp(['Checking for the No. ' num2str(i) ' inversion: Inv ' num2str(pars.Inv(i)) '.'])
            %disp(['This script is for the inversion: No. ' num2str(i) '.'])
            clear scfiles
            Scfiles = {};
            % load all epoch files
            load('/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/EEG/subs-eeg/files.mat')
            % define function
            select_epoch = @(fepoch) (isempty(strfind(fepoch,'green')) & isempty(strfind(fepoch,'red')));
            mask_files   = cell2mat(cellfun(select_epoch, files, 'UniformOutput', false));
            Scfiles      = files(mask_files);
            
            nii_output = {};
            for ss = loopin_s
                sub = ['S' sprintf('%02d', ss)];
                subdir = fullfile(dir_base, sub);
                if w == 1; time_wind = ['-216_1600';'-1516_199'];
                else time_wind = ['-200_*_f'; '-200_*_f']; end
                es_nii = dir(fullfile(subdir, ['es*_' sub '_scan*_' num2str(pars.Inv(i)) '_t' time_wind(1,:) '*_2.nii']));
                eu_nii = dir(fullfile(subdir, ['eu*_' sub '_scan*_' num2str(pars.Inv(i)) '_t' time_wind(2,:) '*_1.nii']));
                es_nii = {es_nii(:).name}';
                eu_nii = {eu_nii(:).name}';
                nii_output(end+1:end+length(es_nii),1) = es_nii;
                nii_output(end+1:end+length(eu_nii),1) = eu_nii;
            end
            
            nii_mat = @(x) [x(1:strfind(x, ['_' num2str(pars.Inv(i)) '_t'])-1) '.mat'];
            done_mat = cellfun(nii_mat, nii_output, 'UniformOutput', false);
            
            select = @(x) x(68:end);
            all_mat = cellfun(select, Scfiles, 'UniformOutput', false);
            
            [torun_files torun_inx] = setdiff(all_mat, done_mat);
            
            if isempty(torun_inx);
                disp(['All results have been written for the freq band: ' num2str(pars.Foi(i,:)) ', time window: ' num2str(Woi(w,:)) '.'])
                continue;
            else
                if min(torun_inx) > 1
                    pudding_file = all_mat(min(torun_inx)-1);
                    torun_files(end+1,1) = pudding_file;
                    torun_inx = [min(torun_inx)-1; torun_inx];
                else
                    torun_inx = torun_inx;
                end
                disp(['Not finished files will be written for the freq band: ' num2str(pars.Foi(i,:)) ', time window: ' num2str(Woi(w,:)) '.'])
                torun_files
                
                %% write matlabbatch
                clear matlabbatch
                matlabbatch{1}.spm.meeg.source.results.D = Scfiles(torun_inx);  % with path, feed to spm matlabbatch
                matlabbatch{1}.spm.meeg.source.results.val = pars.Inv(i); % inversion index
                matlabbatch{1}.spm.meeg.source.results.woi = Woi(w,:);
                matlabbatch{1}.spm.meeg.source.results.foi = pars.Foi(i,:);
                matlabbatch{1}.spm.meeg.source.results.ctype = 'induced'; % single trials
                matlabbatch{1}.spm.meeg.source.results.space = 1;
                matlabbatch{1}.spm.meeg.source.results.format = 'image'; % mni image
                matlabbatch{1}.spm.meeg.source.results.smoothing = 6; % smooth 6 mm
                % run batch
                spm_jobman('run', matlabbatch);
                
            end
        end
    end
end

toc
%diary off
