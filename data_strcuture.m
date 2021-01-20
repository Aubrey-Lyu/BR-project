cd /lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/subs/S01/bhv
load S01_scan1_resp1.mat
load S01_scan1_resp2.mat
load S01_scan1_resp3.mat
load S01_scan1_resp4.mat
load S01_scan1_resp5.mat
load S01_scan1_time_trial.mat
lat_resp = [resp1(:,2)+time_trial(1);...
    resp2(:,2)+time_trial(2);...
    resp3(:,2)+time_trial(3);...
    resp4(:,2)+time_trial(4);...
    resp5(:,2)+time_trial(5)];
code_resp = [resp1(:,1); resp2(:,1); resp3(:,1); resp4(:,1); resp5(:,1)];

% lat_resp2 = [trialResponse(1).domResponse(:,2)+time_trial(1);...
%     trialResponse(2).domResponse(:,2)+time_trial(2);...
%     trialResponse(3).domResponse(:,2)+time_trial(3);...
%     trialResponse(4).domResponse(:,2)+time_trial(4);...
%     trialResponse(5).domResponse(:,2)+time_trial(5)];
% code_resp2 = [trialResponse(1).domResponse(:,1); trialResponse(2).domResponse(:,1); trialResponse(3).domResponse(:,1); trialResponse(4).domResponse(:,1); trialResponse(5).domResponse(:,1)];

cd /lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/subs/S01/
load S01_scan1_filtered_gac_cbc.mat

lat_stim = [latency_s1/250, latency_s2/250, latency_s3/250];
code_stim = [ones(1, length(latency_s1)),2*ones(1, length(latency_s2)),3*ones(1, length(latency_s3))] ;

close all
scatter(lat_resp, code_resp, 'filled', 'MarkerFaceColor',[0 .7 .7])
hold on
scatter(lat_stim, code_stim,...
    'MarkerEdgeColor','r');

hold off
