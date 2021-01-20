clear
load('/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/info/rival_ind.mat');
load('/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/info/replay_smooth_ind.mat'); % index of scans/sessions
load('/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/info/replay_inst_ind.mat');

dir_base = '/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/EEG/subs-eeg';
percepts = {'domin'; 'u'};

%% table variables
Sub = {}; Scan = []; Percept = {}; Cond = {}; Inv = []; F = [];
%%
%spm('defaults', 'EEG');

for ss = 1:20
    sub = ['S' sprintf('%02d', ss)];
    subdir = fullfile(dir_base, sub);
    
    for n = 1:3
        
        if n == 1; IND = rival_ind; cond = 'BR'; 
        elseif n == 2;
            IND = replay_smooth_ind; cond = 'RplSM'; 
        else
            IND = replay_inst_ind; cond = 'RplIN';
        end
        
        for p = 1:length(percepts)
            percept = percepts{p};
            ind = IND{ss};
            for sc = 1: length(ind)
                scfile = fullfile( subdir, ['e' percept '_' sub '_scan' num2str(ind(sc)) '.mat']);
                disp(['Loading for ' scfile]);
                clear D; 
                load(scfile)
                
                for inv = 1:18 
                    F = [F; D.other.inv{1,inv}.inverse.F];
                    Sub{end+1,1} = sub;
                    Scan = [Scan; sc];
                    Percept{end+1,1} = percept;
                    Cond{end+1,1} = cond;
                    Inv = [Inv; inv];
                    disp(['Data has been written for ' 'e' percept '_' sub '_scan' num2str(ind(sc)) '.mat, inversion ' num2str(inv)])
                end
            end 
        end
    end
end

%% write table
T = table();
T.Sub = Sub; T.Scan = Scan; T.Inv = Inv; T.Percept = Percept; T.Cond = Cond; T.F = F;
writetable(T, '/data/dl577/BR/info/F_inv_edominu_allc.csv')