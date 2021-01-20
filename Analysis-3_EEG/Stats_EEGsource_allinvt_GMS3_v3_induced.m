%-----------------------------------------------------------------------
% Source inversion
% modified at 09 Oct 2019 by Dian
% GMS3=> assume independence between conditions within subjects
% no grand mean scalling, but with overall grand mean scaling = 2
% do for all specified inversion models, for all time windows
%-----------------------------------------------------------------
clear
%% SPM set up
rmpath(genpath('/applications/spm/spm12_6906'))
addpath /home/dl577/spm12
addpath ~/scripts/matlab/self-defined_functions/

spm('Defaults', 'fmri');
spm_jobman('initcfg');

%% what analysis
check_nii = 0;
glm1_spec = 0;
glm1_esti = 0;
glm1_cont = 0;
glm2_spec = 0;
glm2_esti = 0;
glm2_nonp_spec = 1;
glm2_nonp_esti = 1;

%% define common key variables - reusable
model_folder = 'model_GMS3_v3';
No_inversion = [2 7 9 11 13 15 17];
Inv = 5 %1:7;
Tw  = 2 %1:2;
loopin = 1:20;
outliers = [];

% loop through all Inv and Tw
for t_index = Tw
    for inv_index = Inv
        % deal with outlier data
        if t_index == 1 && inv_index == 5
            outliers=[13];
        end
        
        inv_model = ['inv' num2str(inv_index) '_t' num2str(t_index)];
        dir_base = ['/data/dl577/BR/EEG_inv_induced/' inv_model];
        
        disp(['Processing for ' inv_model '...'])
        
        glm1_done = 0;
        % reference file = index of good trials
        load('/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/info/good_trials_index_all.mat') % read as the variable "gtIndx"
        mat_lst = gtIndx(:,1);
        Idnx    = gtIndx(:,2);
        
        % Specify contrasts by running script (for common use)
        addpath /home/dl577/scripts/spm_scripts/BR
        run_script_specify_contrasts
        Tvecs = who('tvec*');
        Fvecs = who('fvec*');
        
        % replace special characters in glm2 folder names
        folder_names = strrep(contrasts(1,:),' > ','-');
        folder_names = strrep(folder_names, ' (', '_');
        folder_names = strrep(folder_names, ')', '');
        
        % for checking nifti
        all_sub_not_match_ind = {};
        
        %         % dealing with outliers
        %         if inv_index == 4 && t_index == 1;
        %             loopin = [1:10 12:20];
        %         elseif inv_index == 4 && t_index==2;
        %             loopin = [1:16, 18:20];
        %         else loopin = 1:20;
        %         end
        
        %% Check source image files are complete
        for ss = loopin
            not_match_ind = [];
            sub = ['S' sprintf('%02d', ss)];
            subdir = fullfile(dir_base, sub);
            efn = dir(fullfile(subdir, ['e*scan*_' num2str(No_inversion(inv_index)) '_t*.nii']));
            efn = char({efn(:).name}');
            scfn = unique(cellstr(efn(:, 1:13)));
            
            %% =================== check nii ========================
            good = 0;
            if check_nii
                
                for sc = 1:length(scfn)
                    fn = scfn{sc};
                    if fn(end) == '_'; fn = fn(1:end-1); end
                    
                    disp(['Checking for the scan: ' fn '...'])
                    % find how many nii scans
                    nii_fn = dir(fullfile(subdir, [fn, '_*.nii'])); nii_fn = {nii_fn(:).name}';
                    n_scan = 0;
                    for nn = 1:length(nii_fn)
                        n_scan = n_scan + spm_select_get_nbframes(fullfile(subdir,nii_fn{nn}));
                    end
                    disp(['There are ' num2str(n_scan) ' scans in the result nifti file(s).'])
                    % reference: EEG recorded good trials
                    find_result = strfind(mat_lst, [fn '.mat']);
                    mat_ind = find(cellfun('isempty', find_result)==0);
                    ref_len = length(Idnx{mat_ind});
                    disp(['The original data has ' num2str(ref_len) ' good trials.'])
                    fprintf('\n')
                    
                    % compare and record;
                    % for induced response, there should only be 1 result
                    % across trials
                    
                    if (n_scan ~= 2 && fn(2)=='s') || ((n_scan ~= 1 && fn(2)=='u'))
                        not_match_ind(end+1,1) = mat_ind;
                    end
                end
                all_sub_not_match_ind{ss,1} = not_match_ind;
                
                if isempty(not_match_ind)
                    good = 1;
                    disp('============= Good to go! ==============>>>>>>>')
                    disp(' ')
                else
                    disp('Check your broken files...')
                    save(['/data/dl577/BR/info/all_sub_nii_wrong_induced_' inv_model '.mat'], 'all_sub_not_match_ind')
                    continue
                end
            end
            
            % ----- >  after making sure files are ok, do the following...
            if good || check_nii == 0
                
                %% ============= GLM1 specification ======================
                % --------------------------------------------------------------
                %%%% anova: cell1: BR-green;  cell2: BR-red;  cell3: BR-mixed
                %%%%%%%%%%% cell4: Rpl-green; cell5: Rpl-red; cell6: Rpl-mixed
                %% ====================================================
                if glm1_spec && exist(fullfile(subdir, model_folder, 'SPM.mat'), 'file') ~= 2 % if not yet done
                    disp(['GLM1 specification for the subj: ' sub '.'])
                    output_dir = fullfile(subdir, model_folder);
                    mkdir(output_dir)
                    disp(['Output folder: ' output_dir '.'])
                    % load condition index
                    load('/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/info/rival_ind.mat');
                    load('/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/info/replay_smooth_ind.mat'); % index of scans/sessions
                    load('/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/info/replay_inst_ind.mat');
                    
                    BR_ind = rival_ind{ss}; Rpl_ind = sort([replay_smooth_ind{ss} replay_inst_ind{ss}]);
                    
                    cell1 = {}; cell2 = {}; cell3 = {};
                    cell4 = {}; cell5 = {}; cell6 = {};
                    
                    for sc = 1:length(scfn)
                        fn = scfn{sc};
                        if fn(end) == '_'; fn = fn(1:end-1); end
                        sc_ind = str2num(fn(12:end));
                        
                        %---> for BR cond
                        if ismember(sc_ind, BR_ind)
                            if fn(2) =='s'
                                % for green percept
                                clear nii_fn
                                nii_fn = dir(fullfile(subdir, [fn, '_*1.nii']));
                                nii_3d_fn = cellstr(expand_4d_vols(fullfile(subdir, nii_fn.name)));
                                cell1(end+1:end+length(nii_3d_fn),1) = cellfun(@(x) x(find(~isspace(x))), nii_3d_fn, 'UniformOutput', false);
                                % for red percept
                                clear nii_fn
                                nii_fn = dir(fullfile(subdir, [fn, '_*2.nii']));
                                nii_3d_fn = cellstr(expand_4d_vols(fullfile(subdir, nii_fn.name)));
                                cell2(end+1:end+length(nii_3d_fn),1) = cellfun(@(x) x(find(~isspace(x))), nii_3d_fn, 'UniformOutput', false);
                            elseif fn(2) == 'u'
                                % for mixed percept
                                clear nii_fn
                                nii_fn = dir(fullfile(subdir, [fn, '_*1.nii']));
                                nii_3d_fn = cellstr(expand_4d_vols(fullfile(subdir, nii_fn.name)));
                                cell3(end+1:end+length(nii_3d_fn),1) = cellfun(@(x) x(find(~isspace(x))), nii_3d_fn, 'UniformOutput', false);
                            end
                            %---> for Rpl cond
                        elseif ismember(sc_ind, Rpl_ind)
                            if fn(2) =='s'
                                % for green percept
                                clear nii_fn
                                nii_fn = dir(fullfile(subdir, [fn, '_*1.nii']));
                                nii_3d_fn = cellstr(expand_4d_vols(fullfile(subdir, nii_fn.name)));
                                cell4(end+1:end+length(nii_3d_fn),1) = cellfun(@(x) x(find(~isspace(x))), nii_3d_fn, 'UniformOutput', false);
                                % for red percept
                                clear nii_fn
                                nii_fn = dir(fullfile(subdir, [fn, '_*2.nii']));
                                nii_3d_fn = cellstr(expand_4d_vols(fullfile(subdir, nii_fn.name)));
                                cell5(end+1:end+length(nii_3d_fn),1) = cellfun(@(x) x(find(~isspace(x))), nii_3d_fn, 'UniformOutput', false);
                            elseif fn(2) == 'u'
                                % for mixed percept
                                clear nii_fn
                                nii_fn = dir(fullfile(subdir, [fn, '_*1.nii']));
                                nii_3d_fn = cellstr(expand_4d_vols(fullfile(subdir, nii_fn.name)));
                                cell6(end+1:end+length(nii_3d_fn),1) = cellfun(@(x) x(find(~isspace(x))), nii_3d_fn, 'UniformOutput', false);
                            end % if fn(2)
                        end % if ismember(condition)
                    end % for sc
                    disp('Nifti files have been groups into 6 cells.')
                    disp(['Cell 1-6 have ' num2str(length(cell1)) ', ' num2str(length(cell2)) ', ' num2str(length(cell3)) ', '...
                        num2str(length(cell4)) ', ' num2str(length(cell5)) ', ' num2str(length(cell6)) ' scans, respectively.'])
                    fprintf('\n')
                    if ~ismember(ss, outliers) % if not yet sort out outliers
                        if isempty(cell1) || isempty(cell2) || isempty(cell3) || isempty(cell4) || isempty(cell5) || isempty(cell6)
                            error('One of the cells is empty; deal with the outlier separately.')
                        end
                    end
                    %% write matlabbatch
                    clear matlabbatch
                    matlabbatch{1}.spm.stats.factorial_design.dir = {output_dir};
                    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(1).scans = cell1;
                    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(2).scans = cell2;
                    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(3).scans = cell3;
                    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(4).scans = cell4;
                    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(5).scans = cell5;
                    if ~ismember(ss, outliers) % sub 13 and 11 does not have cell6
                        matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(6).scans = cell6;
                    end
                    matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = 0; % 0: independent; 1: not independent
                    matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = 1; % 1:unequal variance
                    matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0; % global mean scalling
                    matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;
                    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
                    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
                    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
                    matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
                    matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/data/dl577/masks/GM_WM_analysisMask.nii'};
                    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
                    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
                    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 2;
                    % run batch
                    disp('Begin to run...')
                    spm_jobman('run',matlabbatch);
                    disp(['GLM1 specification is done for ' sub '!'])
                end %if glm1
                
                %% ==================== GLM1 estimation ========================
                if glm1_esti && exist(fullfile(subdir, model_folder, 'RPV.nii'), 'file') ~= 2
                    % deal with outlier
                    if ss == 18 && inv_index == 5 
                        continue
                    else
                        spm_file = fullfile(subdir, model_folder,'SPM.mat');
                        clear matlabbatch
                        matlabbatch{1}.spm.stats.fmri_est.spmmat = {spm_file}; %Select the SPM.mat of glm1
                        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
                        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
                        % run batch
                        disp('Begin to run...');
                        spm_jobman('run',matlabbatch);
                        disp(['GLM1 estimation is done for ' sub '!'])
                    end
                end
                
                %% =================   GLM1 contrasts   =======================
                if glm1_cont
                    
                    con_loop = 1: size(contrasts,2);
                    % check if the contrasts are all done (and no more than 15 contrasts)
                    fs_con = dir(fullfile(subdir, model_folder, 'con_0*.nii'));
                    fs_ess = dir(fullfile(subdir, model_folder, 'ess_0*.nii'));
                    %------------------------------------------------------------------%
                    % if the contrast if already completed
                    if length(fs_con) == length(Tvecs) && length(fs_ess) == length(Fvecs);
                        disp(['All GLM1 contrsts are already completed for ' sub ', skipping it...'])
                        continue;
                    else
                        %=================== write matlabbatch ==================
                        clear matlabbatch
                        
                        spm_file = fullfile(subdir, model_folder,'SPM.mat');
                        matlabbatch{1}.spm.stats.con.spmmat = {spm_file};
                        %-------------loop to write matlabbatch parts -------
                        for cc = con_loop
                            cname = contrasts{1,cc};
                            cvec = contrasts{2,cc};
                            %--------------- define contrast --------------
                            if ismember(ss, outliers)
                                if ~isempty(find(cvec(:,6))) % if not applicable contrasts
                                    % create dummy contrasts
                                    cname = 'dummy';
                                    cvec(:) = 0; cvec(1) = 0.01;
                                    % if the 6th colomn is empty, then this
                                    % contrast is applicable to sub13
                                end
                                cvec = cvec(:,1:end-1);
                            end
                            
                            % dealing with t or f contrasts
                            if size(cvec,1) == 1; prefix_con = 'con_'; % if t contrast
                            else prefix_con = 'ess_'; end % or f contrast
                            % for t contrast
                            if size(cvec,1) == 1
                                matlabbatch{1}.spm.stats.con.consess{cc}.tcon.name = cname;
                                matlabbatch{1}.spm.stats.con.consess{cc}.tcon.weights = cvec;
                                matlabbatch{1}.spm.stats.con.consess{cc}.tcon.sessrep = 'none';
                            else % for f contrast
                                matlabbatch{1}.spm.stats.con.consess{cc}.fcon.name = cname;
                                matlabbatch{1}.spm.stats.con.consess{cc}.fcon.weights = cvec;
                                matlabbatch{1}.spm.stats.con.consess{cc}.fcon.sessrep = 'none';
                            end
                        end % end of loop within matlabbatch
                        % continue to write matlabbatch
                        matlabbatch{1}.spm.stats.con.delete = 1; % delete old contrasts to be safe
                        % ---> run batch
                        disp('Begin to run...')
                        spm_jobman('run',matlabbatch);
                        disp([num2str(length(con_loop)) ' contrasts are done for ' sub '!'])
                    end
                    
                end % if glm_cont
            end % if good
        end% for ss
        
        % if completed
        if ss == 20; glm1_done = 1; end
        %---->  GLM2 depends on the completion of GLM1
        if glm1_done == 1 || glm1_cont == 0 % if glm1 is done, do the following...
            
            %% ================== GLM2 specification =======================
            
            if glm2_spec == 1
                
                % deal with outliers
                if  inv_index == 5
                    loopin = [1:17, 19:20];
                else
                    loopin = 1:20;
                end
                
                %================== loop through each contrasts from glm1 =============
                for cc = 1:size(contrasts,2)
                    cname = contrasts{1,cc}; cvec = contrasts{2,cc};
                    glm2_folder = fullfile(dir_base, model_folder, [folder_names{cc} '_' num2str(cc)]);
                    
                    % dealing with t or f contrasts
                    if size(cvec,1) == 1; prefix_con = 'con_'; % if t contrast
                    else prefix_con = 'ess_'; end % or f contrast
                    
                    % now the "sub" is the last subject from the previous glm1
                    if exist(fullfile(glm2_folder,'SPM.mat'),'file') ~= 2 % if these files do not exist
                        
                        % dealing with outlier sub
                        if isempty(find(cvec(:,6))); loopin = loopin;
                        else
                            [~,I_out,~]=intersect(loopin, outliers);
                            loopin(I_out) = [];
                        end
                        
                        %------------------ write matlabbatch ----------------------------------
                        clear matlabbatch
                        
                        mkdir(glm2_folder);
                        matlabbatch{1}.spm.stats.factorial_design.dir = {glm2_folder};
                        % prepare scans
                        scans = {};
                        for ss = loopin
                            sub      = ['S' sprintf('%02d',ss)];
                            GLM1_dir = fullfile(dir_base, sub, model_folder);
                            scans{end+1,1} = fullfile(GLM1_dir, [prefix_con sprintf('%04d', cc) '.nii,1']);
                        end
                        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = scans;
                        % other parameters
                        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
                        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
                        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
                        matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
                        matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/data/dl577/masks/GM_WM_analysisMask.nii'};
                        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
                        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
                        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
                        % run
                        spm_jobman('run',matlabbatch);
                        disp(['One sample T done for the contrast: ' cname '.'])
                    end % if glm1 contrast's been done for the last sub
                end
            end
            
            %% ===================== glm2 estimation ===========================
            if glm2_esti == 1

                for cc = 1:size(contrasts,2)
                    cvec = contrasts{2,cc};
                    glm2_folder = fullfile(dir_base,model_folder, [folder_names{cc} '_' num2str(cc)]);
                    if exist(fullfile(glm2_folder, 'ResMS.nii'), 'file') ~= 2
                        clear matlabbatch
                        matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(glm2_folder,'SPM.mat')}; %Select the SPM.mat of glm2
                        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
                        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
                        % run
                        spm_jobman('run',matlabbatch);
                        disp(['One sample T esitimated for the contrast: ' cname '.'])
                    end % if esti's not yet done
                end
            end
            
            %% =============== glm2 non-parametrical testing =======================
            if glm2_nonp_spec
                                
                % deal with outliers
                if  inv_index == 5
                    loopin = [1:17, 19:20];
                else
                    loopin = 1:20;
                end
                
                %-----------loop through each contrasts from glm1 -----
                for cc = 1:size(contrasts,2)
                    cname = contrasts{1,cc}; cvec = contrasts{2,cc};
                    glm2np_folder = fullfile(dir_base, [model_folder '_NP'], [folder_names{cc} '_' num2str(cc)]);
                    % dealing with t or f contrasts
                    if size(cvec,1) == 1; prefix_con = 'con_'; % if t contrast
                    else prefix_con = 'ess_'; end % or f contrast
                    
                    % now the "sub" is the last subject from the previous glm1
                    if exist(fullfile(glm2np_folder,'SnPMcfg.mat')) ~= 2 % while glm specification is not done
                        
                        % dealing with outlier sub
                        if isempty(find(cvec(:,6))); loopin = loopin;
                        else
                            [~,I_out,~]=intersect(loopin, outliers);
                            loopin(I_out) = [];
                        end
                        %------------------ write matlabbatch ----------------------------------
                        clear matlabbatch
                        
                        mkdir(glm2np_folder);
                        
                        % prepare scans
                        scans = {};
                        for ss = loopin
                            sub      = ['S' sprintf('%02d',ss)];
                            GLM1_dir = fullfile(dir_base, sub, model_folder);
                            scans{end+1,1} = fullfile(GLM1_dir, [prefix_con sprintf('%04d', cc) '.nii,1']);
                        end
                        matlabbatch{1}.spm.tools.snpm.des.OneSampT.DesignName = 'MultiSub: One Sample T test on diffs/contrasts';
                        matlabbatch{1}.spm.tools.snpm.des.OneSampT.DesignFile = 'snpm_bch_ui_OneSampT';
                        matlabbatch{1}.spm.tools.snpm.des.OneSampT.dir = {glm2np_folder};
                        matlabbatch{1}.spm.tools.snpm.des.OneSampT.P = scans;
                        matlabbatch{1}.spm.tools.snpm.des.OneSampT.cov = struct('c', {}, 'cname', {});
                        matlabbatch{1}.spm.tools.snpm.des.OneSampT.nPerm = 5000;
                        matlabbatch{1}.spm.tools.snpm.des.OneSampT.vFWHM = [0 0 0];
                        matlabbatch{1}.spm.tools.snpm.des.OneSampT.bVolm = 1;
                        matlabbatch{1}.spm.tools.snpm.des.OneSampT.ST.ST_U = 2.03;
                        matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.tm.tm_none = 1;
                        matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.im = 0;
                        matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.em = {'/data/dl577/masks/GM_WM_analysisMask.nii'};;
                        matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalc.g_omit = 1;
                        matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.gmsca.gmsca_no = 1;
                        matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.glonorm = 2;
                        % run
                        spm_jobman('run',matlabbatch);
                        disp(['Non-parametric one sample T specified for the No.' num2str(cc) ' contrast: ' cname '.'])
                    end % if exist
                end % for cc
            end
            
            %% compute glm2 non-parametric test
            if glm2_nonp_esti
                for cc = 1:size(contrasts,2)
                    cname = contrasts{1,cc}; cvec = contrasts{2,cc};
                    glm2np_folder = fullfile(dir_base, [model_folder '_NP'], [folder_names{cc} '_' num2str(cc)]);
                    
                    if exist(fullfile(glm2np_folder,'ResMS.img')) ~= 2
                        %------------------ write matlabbatch ----------------------------------
                        clear matlabbatch
                        matlabbatch{1}.spm.tools.snpm.cp.snpmcfg = {fullfile(glm2np_folder, 'SnPMcfg.mat')};
                        % run
                        spm_jobman('run',matlabbatch);
                        disp(['Non-parametric one sample T esitimated for the No.' num2str(cc) ' contrast: ' cname '.'])
                    end
                end
            end % if to compute non-p glm2
            %--------------------------------------------------------------------------------------------
        else % if GLM1 is done
            continue % if not, continue to the next inversion model. (save time when running a script)
        end % if GLM1 is done
    end % for inv index
end % for tw