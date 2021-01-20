%-----------------------------------------------------------------------
% Job saved on 01-Aug-2019 20:55:07 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (6906)
% cfg_basicio BasicIO - Unknown
% I crewed up the trial events (and consequentially the labelling of bad
% epochs) by mistakenly adding a padding, so there is an additional
% procedure of "corr_events".
% correct trial events are stored in the path: /lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/EEG/SPM_objects
%-----------------------------------------------------------------------
clear
rmpath(genpath('/applications/spm/spm12_6906'))
addpath('/home/dl577/spm12')
addpath ~/scripts/matlab/self-defined_functions/
spm('defaults', 'EEG');
spm_jobman('initcfg');

dir_base = '/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/EEG/subs-eeg';
cd(dir_base)


%% what analyses?
epoch = 1;
rdgn  = 1;
red   = 0;
green = 0;
mixed = 0;
clean = 0;
corr_events = 0; % after messing up with padding...

%%

if epoch
    sfiles = {};
%     for ss = 1:20
%         sub = ['S' sprintf('%02d%',ss)];
%         
%     sfs = dir(fullfile(dir_base, sub, [sub '_scan*.mat']));
%     for i = 1:length(sfs)
%         sfiles{end+1,1} = fullfile(dir_base, sub, sfs(i).name);
%     end
%     end
   sfiles = {fullfile(dir_base, 'S10', 'S10_scan7.mat');
        
        };
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:size(sfiles,1)
        file = sfiles{i};
        load(file); % read as D
        D.channels(32).bad = 1; % ECG
        %     pnt = D.fiducials.fid.pnt;
        %     pnt_new = [pnt(:,2) pnt(:,1) pnt(:,3)];
        %     D.fiducials.fid.pnt = pnt/10;
        %     D.sensors.eeg.chanpos = pnt/10;
        %     D.sensors.eeg.elecpos = pnt/10;
        %     D.fiducials.fid.pnt = pnt_new;
        %     D.sensors.eeg.chanpos = pnt_new;
        %     D.sensors.eeg.elecpos = pnt_new;
        %
        save(file, 'D')
        
        fn = strsplit(file,'/');
        sub = fn{end-1};
        matfile = fn{end};
        subdir = fullfile(dir_base, sub);
        cd(subdir) % subj folder
        
        %% for stable percepts
        % in the version2, change the epoch time
        if rdgn == 1
            clear matlabbatch
            matlabbatch{1}.spm.meeg.preproc.epoch.D = {file};
            matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.timewin = [-1 0.1]*1000; % in ms
            matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).conditionlabel = 'stable_green';
            matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventtype = 'Stimulus';
            matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventvalue = 'S  1';
            matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).trlshift = 0;
            matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).conditionlabel = 'stable_red';
            matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).eventtype = 'Stimulus';
            matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).eventvalue = 'S  2';
            matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).trlshift = 0;
            matlabbatch{1}.spm.meeg.preproc.epoch.bc = 0;
            matlabbatch{1}.spm.meeg.preproc.epoch.eventpadding = 0;%-0.5;
            matlabbatch{1}.spm.meeg.preproc.epoch.prefix = 'edomin_';
            
            spm_jobman('run', matlabbatch);
            disp(['Epoch done for stable percepts for ' matfile])
        end
        
        %% for red percepts
        if red == 1
            clear matlabbatch
            matlabbatch{1}.spm.meeg.preproc.epoch.D = {file};
            matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.timewin = [-1 0.1]*1000; % in ms
            matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).conditionlabel = 'stable_red';
            matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventtype = 'Stimulus';
            matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventvalue = 'S  2';
            matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).trlshift = 0;
            
            matlabbatch{1}.spm.meeg.preproc.epoch.bc = 0;
            matlabbatch{1}.spm.meeg.preproc.epoch.eventpadding = 0;%-0.5;
            matlabbatch{1}.spm.meeg.preproc.epoch.prefix = 'ered_';
            
            spm_jobman('run', matlabbatch);
            disp(['Epoch done for red percepts for ' matfile])
        end
    end
    %% for green percepts
    if green == 1
        clear matlabbatch
        matlabbatch{1}.spm.meeg.preproc.epoch.D = {file};
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.timewin = [-1 0.1]*1000; % in ms
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).conditionlabel = 'stable_green';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventtype = 'Stimulus';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventvalue = 'S  1';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).trlshift = 0;
        
        matlabbatch{1}.spm.meeg.preproc.epoch.bc = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.eventpadding = 0;%-0.5;
        matlabbatch{1}.spm.meeg.preproc.epoch.prefix = 'egreen_';
        
        spm_jobman('run', matlabbatch);
        disp(['Epoch done for green percepts for ' matfile])
    end
    
    
    %% for unstable percepts
    if mixed == 1
        clear matlabbatch
        matlabbatch{1}.spm.meeg.preproc.epoch.D = {file};
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.timewin = [-1 0.1]*1000;
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).conditionlabel = 'unstable_mixed';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventtype = 'Stimulus';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventvalue = 'S  3';
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).trlshift = 0;
        
        matlabbatch{1}.spm.meeg.preproc.epoch.bc = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.eventpadding = 0;%0.5;
        matlabbatch{1}.spm.meeg.preproc.epoch.prefix = 'eu_';
        
        spm_jobman('run', matlabbatch);
        disp(['Epoch done for unstable percept for ' matfile])
    end
end % end of epoch

%% correct trial events
if corr_events == 1
    load('/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/EEG/subs-eeg/files.mat')
    
    for f = 1:length(files)
        clear D
        
        file = files{f};
        fns  = strsplit(file, '/'); efn = fns{end};
        subdir = strjoin({fns{1:end-1}},'/');
        
        cd(subdir); load(file)
        
        clear D_ref
        D_ref  = load(['/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/EEG/SPM_objects/' efn]);
        events = {D_ref.D.trials.events};
        
        [D.trials.events] = events{:};
        
        save(file, 'D')
        disp(['Epoch events corrected for ' efn '.'])
    end
end

%% clean epoch (assign bad to those with more than single event)
if clean == 1
   % load('/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/EEG/subs-eeg/files.mat')
    %files = files(271);
%     % modified at 2 Sep, for rewriting broken files:
%     broken_subs = {'S06', 'S11'};
% files = {};
% for bs = 1:length(broken_subs)
%     broken_sub = broken_subs{bs};
%     broken_files = dir(fullfile(dir_base, broken_sub, 'es*.mat'));
%     for bf = 1:length(broken_files)
%         broken_file = broken_files(bf).name;
%         files{end+1,1} = fullfile(dir_base, broken_sub, broken_file);
%     end
% end
% files = {fullfile(dir_base, 'S06', 'es_S06_scan3.mat');
%         fullfile(dir_base, 'S11', 'es_S11_scan11.mat');
%         fullfile(dir_base, 'S06', 'es_S06_scan10.mat');
%         };
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    files = {fullfile(dir_base, 'S01', 'es_S01_scan1.mat')};
    for f = 1:length(files)
        clear D
        
        file = files{f};
        fns  = strsplit(file, '/'); efn = fns{end};
        subdir = strjoin({fns{1:end-1}},'/');
        
        cd(subdir); load(file)
        
        clear matlabbatch D
        
        load(file)
        
        for i = 1:length(D.trials)
            D.trials(i).bad = 0; % default is 0
            if ~isempty(D.trials(i).events)
                if length(cell2mat(strfind({D.trials(i).events.code}, 'Stimulus'))) > 1 % having more stimuli than the onset
                    D.trials(i).bad = 1;
                end
            end
        end
        
        save(file, 'D')
        disp(['Epoch cleaned for ' efn '.'])
    end
end

