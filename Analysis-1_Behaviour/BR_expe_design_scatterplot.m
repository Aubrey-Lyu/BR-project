%% model individual behavioural responses
clear
addpath /data/dl577/scripts/matlab/
close all
ax = axes;
coltab =  [11,119,94; 244,48,15; 210,210,210 ]; % green, red, white ;
%yellow: [253, 210,98]
%coltab =  [70,172,200; 180,15,32 ; 210,210,210 ]; % blue, red, white
altab  = [1; 0.3 ;1];
scans_info = [13,12,12,12,12,13,12,9,9,12,13,17,10,12,12,11,14,12,12,12];
xlim([0, 32])
ylim([0,280])

n = 0; % counter
for ss = 1:20
    subdir = ['/lustre/scratch/wbic-beta/dl577/Binocular_Rivalry/subs/S' num2str(ss, '%02.f') '/bhv'];
    brk = 25/sum(scans_info);%1/scans_info(ss);
    n = n+0.7*brk;
    vline(n , 'y' ,['S' num2str(ss, '%02.f')]);
    for sc = 1:scans_info(ss)
        load([subdir '/S' num2str(ss, '%02.f') '_scan' num2str(sc) '_time_trial.mat'])
        load([subdir '/S' num2str(ss, '%02.f') '_scan' num2str(sc) '_resp1.mat'])
        load([subdir '/S' num2str(ss, '%02.f') '_scan' num2str(sc) '_resp2.mat'])
        load([subdir '/S' num2str(ss, '%02.f') '_scan' num2str(sc) '_resp3.mat'])
        load([subdir '/S' num2str(ss, '%02.f') '_scan' num2str(sc) '_resp4.mat'])
        load([subdir '/S' num2str(ss, '%02.f') '_scan' num2str(sc) '_resp5.mat'])
        rp1_choise = resp1(:,1);
        rp1_time = resp1(:,2) + time_trial(1);
        rp2_choise = resp2(:,1);
        rp2_time = resp2(:,2) + time_trial(2);
        rp3_choise = resp3(:,1);
        rp3_time = resp3(:,2) + time_trial(3);
        rp4_choise = resp4(:,1);
        rp4_time = resp4(:,2) + time_trial(4);
        rp5_choise = resp5(:,1);
        rp5_time = resp5(:,2) + time_trial(5);
        
        c = [rp1_choise; rp2_choise;rp3_choise;rp4_choise;rp5_choise];
        
        RGB = coltab(c(:), :)/256;
     %   ALPHA = altab(c(:), :);
        n = n + brk;
        t = [rp1_time; rp2_time;rp3_time;rp4_time;rp5_time];
        x = ones(1, length(t))*n;
        
        %% make plots
        
        hold on;
       % for i = 1:length(t)
            s = scatter(x,t, 4, RGB, 's'); set(gca,'Color','k')
            s.MarkerFaceAlpha = 1;
          %  s.MarkerFaceAlpha = ALPHA(i);
        %end
        %
        % plot(rp1_time, rp1_choise)
        % hold on
        % plot(rp2_time, rp2_choise)
        % hold on
        % plot(rp3_time, rp3_choise)
        % hold on
        % plot(rp4_time, rp4_choise)
        % hold on
        % plot(rp5_time, rp5_choise)
    end
   
    
end
