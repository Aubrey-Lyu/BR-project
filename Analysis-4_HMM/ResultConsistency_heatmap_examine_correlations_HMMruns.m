%% plot heatmap for HMM dignosis of result stability

result_folder = '/media/dian/D/data/Binocular_Rivalry/HMM_analyses/output_test/embedded_lag/DMN_V1_1000_100_K4_pca16/run1_analyses';
cd(result_folder)
load('examine_state_correlation.mat')



%% relabel column names
labels = cell(length(colname),1);
for i = 1:length(colname)
    c = colname{i};
    c = strrep(c, 'half','Half');
    c = strrep(c, 'full','Full ');
    keys = strsplit(c, '_');
    redo_c = [keys{2}, ' ', keys{1}];
    labels{i} = redo_c;
end


%% heat map
close all
figure;heatmap(correlation_MAT, 'Colormap', summer, 'XData', labels, 'YData', labels)


%% stats
R_diff = zeros(size(v_meanR_not_match));
for i = 1:size(v_Rmax, 2)
    R_diff(:,i) = v_meanR_not_match(:,i) - v_Rmax(1,i);
end

r_mean_diff = mean(R_diff,2);

%% histgram

