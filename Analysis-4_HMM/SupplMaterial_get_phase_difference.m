load /data/dl577/BR/HMM_analyses/embedded_lag/DMN_V1_1000_100_K4_pca16/run1_analyses/toShruti/data.mat
% continuum requencies, stats with all subjects grouped together.
% with the help of the online code: https://uk.mathworks.com/matlabcentral/answers/91647-how-do-i-calculate-the-amplitude-ratio-and-phase-lag-for-two-sinusoidal-signals-in-matlab

ttrial = T(1);
conditions = {'BR_dominant', 'BR_mixed', 'Rpl_smt_dominant', 'Rpl_smt_mixed'};
VOIs = {'rHP', 'lHP', 'lIPL', 'rIPL', 'PCU','PCC', 'ACC', 'V1'};

% %---- generate random color -------
% % Define the two colormaps.
% cmap1 = hot(8);
% cmap2 = winter(8) ;
% % Combine them into one tall colormap.
% combinedColorMap = [cmap1; cmap2];
% % Pick 15 rows at random.
% randomRows = randi(size(combinedColorMap, 1), [8, 1]);
% % Extract the rows from the combined color map.
% randomColors = combinedColorMap(randomRows, :);
% % ---------------------------------------
randomColors = [0, 74, 61; 31, 205, 109; 43, 150, 222; 156, 85, 183;...
    247, 219, 105; 231, 125, 4;193, 55, 35]/256;

% Sampling frequency
Fs = 250;

x = reshape(X, ttrial, size(X,1)/ttrial, size(X,2));
y = reshape(Y, ttrial, size(X,1)/ttrial, size(Y,2));

close all

for c = 1:4
    condition = conditions{c};
    figure;
    cI = find(y(1,:, c+4)==1);
    
    cy = y(:,cI,:);
    cx = x(:,cI,:);
    
    N_s = length(cI);
    PHASE_LAG = zeros(N_s,40);
    
    for v = 1:length(VOIs)-1
        for i = 1:N_s
            s0 = squeeze(cx(:,i,8)); % v1
            s1 = squeeze(cx(:,i,v));
            
            % center
            s0 = s0 - mean(s0);
            s1 = s1 - mean(s1);
            
            % take the FFT
            sf0=fft(s0);
            sf1=fft(s1);
            % Calculate the numberof unique points
            NumUniquePts = ceil((ttrial+1)/2);
            %
            % figure
            % subplot(211);
            % f = (0:NumUniquePts-1)*Fs/ttrial;
            % plot(f,abs(sf0(1:NumUniquePts)));
            % title('sf0(f) : Magnitude response');
            % ylabel('|sf0(f)|')
            % subplot(212)
            % plot(f,abs(sf1(1:NumUniquePts)));
            % title('sf1(f) : Magnitude response')
            % xlabel('Frequency (Hz)');
            % ylabel('|sf1(f)|')
            % figure(3)
            % subplot(211)
            % plot(f,angle(sf0(1:NumUniquePts)));
            % title('sf0(f) : Phase response');
            % ylabel('Phase (rad)');
            % subplot(212)
            % plot(f,angle(sf1(1:NumUniquePts)));
            % title('sf1(f) : Phase response');
            % xlabel('Frequency (Hz)');
            % ylabel('Phase (rad)');
            % Determine the max value and max point.
            % This is where the sinusoidal
            %             % is located. See Figure 2.
            %             [mag_x idx_x] = max(abs(sf0));
            %             [mag_y idx_y] = max(abs(sf1));
            %
            %             % determine the phase difference
            %             % at the maximum point.
            %             px = angle(sf0(idx_x));
            %             py = angle(sf1(idx_y));
            %             phase_lag = py - px;
            
            % determine the phase difference
            % at each frequency.
            px = angle(sf0(1:NumUniquePts));
            py = angle(sf1(1:NumUniquePts));
            phase_lag = py - px;
            
            PHASE_LAG(i,:) = phase_lag(1:40)';
            
           
        end % for i = 1:N_s trial
         h = ttest(PHASE_LAG);
            mean_PHASE_LAG = mean(PHASE_LAG);
            % frequency with significant (stable) phase coherence
            fI = find(h==1);
            sig_PHASE_LAG = mean_PHASE_LAG(fI);
            
            scatter(fI, sig_PHASE_LAG, 'MarkerEdgeColor', randomColors(v,:),...
                'MarkerFaceColor', randomColors(v,:), 'MarkerFaceAlpha',.5); hold on
    end % for region
    
    
    hold off
%     set(gca,'YTick', -pi:pi/2:pi)
%     set(gca,'YTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
    legend(VOIs{1:end-1})
end % for condition



