s1 ={'rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth','replay_inst','rivalry','replay_smooth','rivalry','replay_inst'};

s2 = {'rivalry','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_inst','rivalry','replay_smooth','replay_inst ','rivalry','replay_smooth'};

s3 = {'rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth','rivalry',...
    'replay_inst','replay_smooth','replay_inst'};

s4 = {'rivalry','replay_smooth','rivalry','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_inst',...
    'replay_smooth','replay_inst','replay_smooth'};

s5 = {'rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth',...
    'rivalry','replay_smooth','replay_inst','rivalry','replay_smooth','replay_inst'};

s6 = {'rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_inst','replay_smooth',...
    'rivalry','replay_smooth','replay_inst','rivalry'};

s7 = {'rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth','rivalry',...
    'replay_inst','replay_smooth','replay_inst'};

s8 = {'rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_inst','rivalry'};

s9 = {'rivalry','rivalry','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_inst','replay_inst'};

s10 = {'rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth',...
    'rivalry','replay_smooth','replay_inst','rivalry','replay_inst','replay_smooth'};

s11 = {'rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth','rivalry',...
    'rivalry','replay_smooth','replay_inst','rivalry','replay_inst','replay_smooth'};

s12 = {'rivalry','replay_smooth','rivalry','rivalry','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth',...
'rivalry','replay_smooth','replay_inst','rivalry','replay_smooth','replay_inst','rivalry'};

s13 ={'rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth','rivalry',...
    'replay_smooth','replay_inst','rivalry'};

s14 = {'rivalry','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_inst',...
    'replay_smooth','replay_smooth','replay_inst'};

s15 = {'rivalry','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth',...
    'replay_inst','replay_smooth','replay_inst'};

s16 = {'replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth',...
    'rivalry','replay_inst'};

s17 = {'rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth','rivalry',...
    'replay_smooth','replay_inst','rivalry','replay_smooth','replay_inst'};

s18 = {'rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_inst','replay_smooth',...
    'rivalry','replay_smooth','replay_inst'};

s19 = {'rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_inst','rivalry',...
    'replay_smooth','replay_inst','replay_smooth'};

s20 = {'rivalry','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth','rivalry','replay_smooth',...
    'replay_inst','replay_smooth','replay_inst'};

SC = {s1, s2, s3, s4, s5, s6, s7, s8, s9, s10,...
    s11, s12, s13, s14, s15, s16, s17, s18, s19, s20};



%% ===================== record index ========================
% wrapper = @(x) find(strcmp(x, 'rivalry'));
% 
% rival_ind = cellfun(wrapper, SC,'UniformOutput', false );
% 
% wrapper = @(x) find(strcmp(x, 'replay_smooth'));
% replay_smooth_ind = cellfun(wrapper, SC,'UniformOutput', false ); 
% 
% wrapper = @(x) find(strcmp(x, 'replay_inst'));
% replay_inst_ind = cellfun(wrapper, SC,'UniformOutput', false );
% %save('/media/dian/D/R/project/BR/data/rivalind.mat', 'rivalind')
% 
% save('/media/dian/D/R/project/BR/data/info/rival_ind.mat', 'rival_ind')
% save('/media/dian/D/R/project/BR/data/info/replay_smooth_ind.mat', 'replay_smooth_ind')
% save('/media/dian/D/R/project/BR/data/info/replay_inst_ind.mat', 'replay_inst_ind')




%% =====================just record boolin========================
% wrapper = @(x) strcmp(x, 'rivalry');
% 
% rival_ind = cellfun(wrapper, SC,'UniformOutput', false );
% 
% wrapper = @(x) strcmp(x, 'replay_smooth');
% replay_smooth_ind = cellfun(wrapper, SC,'UniformOutput', false );
% 
% wrapper = @(x) strcmp(x, 'replay_inst');
% replay_inst_ind = cellfun(wrapper, SC,'UniformOutput', false );
% %save('/media/dian/D/R/project/BR/data/rivalind.mat', 'rivalind')
% 
% save('/media/dian/D/R/project/BR/data/info/rival_bool.mat', 'rival_ind')
% save('/media/dian/D/R/project/BR/data/info/replay_smooth_bool.mat', 'replay_smooth_ind')
% save('/media/dian/D/R/project/BR/data/info/replay_inst_bool.mat', 'replay_inst_ind')
