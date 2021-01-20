function ICN_atl_visualise(positives,negatives,type_of_chart,color_mode,graph_labels,graph_title,ecc_ticks,tick_orientation)
% ICN_atl_visualise(positives,negatives,type_of_chart,color_mode)
% function to visualise some (preferably ICN_atlas) data in rose/polar/bar plots
%
% ICN_atl_visualise is part of the ICN_atlas toolbox
%
% Please use ICN_atlas.m for full functionality
%
% =========================================================================
% THIS IS STILL WORK IN PROGRESS
% =========================================================================
%
% positives: a 1xn vector or a 2xn matrix or and empty vector
%   this data is presented in reddish/yellowish/lighter colors or bright grays
%   1xn vector: values presented as distance from 0
%   2xn matrix: first row codes distance from 0, secon row codes error
%       bars if the latter are shown on a given type of chart
%   [] vector: no positively marked/colored values are shown
%
% negatives: a 1xn vector or a 2xn matrix or and empty vector
%   this data is presented in blueish/greenish/darker colors or dark grays
%   1xn vector: values presented as distance from 0
%   2xn matrix: first row codes distance from 0, secon row codes error bars
%       if the latter are shown on a given type of chart
%   [] vector: no negatively marked/colored values are shown
%
% type_of_chart:
%   1: rose (wedge) plot with centerpoint symmetry between positives and negatives
%   2: rose (wedge) plot with left-right symmetry between positives and negatives
%   3: line plot with centerpoint symmetry between positives and negatives
%   4: line plot with left-right symmetry between positives and negatives
%   5: polar (patch) plot with centerpoint symmetry between positives and negatives
%   6: polar (patch) plot with left-right symmetry between positives and negatives
%   7: polar (patch) plot with overlap between positives and negatives
%   8: filled polar (patch) plot with centerpoint symmetry between positives and negatives
%   9: filled polar (patch) plot with left-right symmetry between positives and negatives
%  10: vertical bars
%  11: horizontal bars
%  12: positives: filled rose (wedge), negatives unfilled rose (wedge) with overlap
%
% color_mode:
%   1: red-blue
%   2: grayscale
%   3: proportional (jet)
%   4: proportional  (parula)
%   5: proportional (KLR1) % colormap(vertcat([0+(1:-0.01:0)'/2 0+(1:-0.01:0)'/1 ones(size((0:0.01:1)'))],[ones(size((0:0.01:1)')) 0+(0:0.01:1)'/1 0+(0:0.01:1)'/3]))
%   6: proportional (KLR2) % colormap(vertcat([0+(1:-0.01:0)'/3 0+(1:-0.01:0)'/2 0.5*ones(size((1:-0.01:0)'/2))],[0.5*ones(size((0:0.01:1)'))/2 0+(0:0.01:1)'/2 0+(0:0.01:1)'/3]))
%   7: cognitive domain-based
%
%__________________________________________________________________________
% 2012-18 Lajos R Kozak MD, PhD / icnatlas.toolbox@gmail.com
% Semmelweis University MR Research Center, Budapest, Hungary
%
% ICN_atl_visualise.m is part of the ICN_Atlas toolbox 
%
% ICN_Atlas is an analysis toolbox aimed at facilitating the interpretation
% of fMRI data in the context of intrinsic connectivity networks (ICN). It 
% was developed by members of the MR Research Center, Semmelweis University,
% Budapest, Hungary and the Institute of Neurology, University College 
% London, UK. Original idea: Louis Lemieux, software development: Lajos 
% R Kozak, testing & validation: Louis Andre van Graan, Adam Szabo, Umair J
% Chaudhary and others.
%
% The ICN_Atlas toolbox is (C) by Lajos R. Kozak & Louis Lemieux under a 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% License (CC-BY-NC-SA-4.0): http://creativecommons.org/licenses/by-nc-sa/4.0/
% For uses or applications not permitted by the CC-BY-NC-SA-4.0 license 
% please contact the authors via e-mail at icnatlas.toolbox@gmail.com 
%
% When using the toolbox please cite: LR Kozak, LA van Graan, UJ Chaudhary, 
% AG Szabo, L Lemieux: ICN_Atlas: Automated description and quantification 
% of functional MRI activation patterns in the framework of intrinsic 
% connectivity networks, NeuroImage 163 (2017) 319-341, 
% DOI: https://doi.org/10.1016/j.neuroimage.2017.09.014
%
% Please see the documents in the LICENSES subdirectory for (a) the license
% terms (ICN_Atlas-license-CC-BY-NC-SA-4.0-legalcode.txt) and (b) for the 
% copyright notices (ICN_Atlas-copyright-notices.txt) detailing information
% on copyright and attribution related to ICN_Atlas, and (c) for the 
% licenses and/or copyright notices related to other software tools, 
% atlases, etc. that were used during the development of this toolbox 
%
% You should have received a copy of the license along with this work, as 
% outlined above. If not, see http://creativecommons.org/licenses/by-nc-sa/4.0/
%
% The full public release version(s) of the unmodified ICN_Atlas toolbox 
% are available at https://www.nitrc.org/projects/icn_atlas/ and/or 
% http://icnatlas.com . For development/pre-release versions and possible 
% collaborations please contact the authors at icnatlas.toolbox@gmail.com
%

% $Id: ICN_atl_visualise.m 86 2018-03-06 15:07:09Z LRKozak $
%
% =========================================================================
% THIS IS STILL WORK IN PROGRESS, please report all bugs encountered
% =========================================================================

% default line width for the polar plots
polarlinewidth=0.5;

% text orientation in polar plots: 0 = parallel, 1 = perpendicular to radii
if nargin<8
    polar_text_perpendicular=0;
else
    polar_text_perpendicular=tick_orientation;
end

% defaults ticks 
if nargin<7
    ecc_ticks=[0.25 0.5 0.75 1 1];
end

% tick can be defined in various ways
switch length(ecc_ticks)
    case 1 % sigle value given 
        if ecc_ticks>0
            % graph limit is at the defined single value
            forced_max=ecc_ticks;
        else
            % graph limit is at the max of data
            forced_max=0;
        end
        % ticks at 25%, 50%, and 75% of max
        forced_ticks=[0.25 0.5 0.75 1];
    case 2 % a two-element vector is given
        if ecc_ticks(2)>0
            % graph limit is at the value of the 2nd element
            forced_max=ecc_ticks(2);
        else
            % if the 2nd element is 0 then the graph limit is at the max of data
            forced_max=0;
        end
        % number of ticks given in the first element of the vector
        % forced ticks are in the 0:1 range
        forced_ticks=0:(1/(ecc_ticks(1)+1)):1;
        forced_ticks=forced_ticks(2:(end));
    otherwise
        % multi-element vector
        % max is at the value of the last element
        forced_max=ecc_ticks(end);
        if forced_max<ecc_ticks(end-1)
            forced_max=ecc_ticks(end-1);
        end
        % ticks are at the defined values (they are projected to the 0:1
        % interval below)
        forced_ticks=ecc_ticks(1:end)/ecc_ticks(end);
        forced_ticks=forced_ticks(1:(end-1));
        ft_OK=find(forced_ticks>0);
        forced_ticks=forced_ticks(ft_OK);
        forced_ticks=sort(forced_ticks);
end

% default graph title is empty
if nargin<6
    graph_title='';
end

% if nargin<5 ==> see below, at %% arrange labels

% default color mode is red-blue
if nargin<4
    color_mode=1;%red-blue
    %color_mode=2;%grayscale
    %color_mode=3;%proportional  jet
    %color_mode=4;%proportional  parula
    %color_mode=5;%proportional KLR1
    %color_mode=6;%proportional KLR2
    %color_mode=7;%cognitive domain-based
end

% default chart is rose (wedge) plot with centerpoint symmetry between positives and negatives
if nargin<3
    type_of_chart=1;
    % type_of_chart:
    %   1: rose (wedge) plot with centerpoint symmetry between positives and negatives
    %   2: rose (wedge) plot with left-right symmetry between positives and negatives
    %   3: line plot with centerpoint symmetry between positives and negatives
    %   4: line plot with left-right symmetry between positives and negatives
    %   5: polar (patch) plot with centerpoint symmetry between positives and negatives
    %   6: polar (patch) plot with left-right symmetry between positives and negatives
    %   7: polar (patch) plot with overlap between positives and negatives
    %   8: filled polar (patch) plot with centerpoint symmetry between positives and negatives
    %   9: filled polar (patch) plot with left-right symmetry between positives and negatives
    %  10: vertical bars
    %  11: horizontal bars
    %  12: positives: filled rose (wedge), negatives unfilled rose (wedge) with overlapend
end

if nargin<2
    % example data for testing purposes
    negatives=[  0.8147    0.9058    0.1270    0.9134    0.6324    0.0975    0.2785    0.5469    0.9575    0.9649];
    negatives(2,:)=[  0.8147    0.9058    0.1270    0.9134    0.6324    0.0975    0.2785    0.5469    0.9575    0.9649]/10;
end
if nargin<1
    % example data for testing purposes
    positives=[ 0.1576    0.9706    0.9572    0.4854    0.8003    0.1419    0.4218    0.9157    0.7922    0.9595];
    positives(2,:)=[ 0.1576    0.9706    0.9572    0.4854    0.8003    0.1419    0.4218    0.9157    0.7922    0.9595]/10;
end
if isempty(negatives) && isempty(positives)
    disp('No data to visualise, exiting.')
    return
end
%% arrange labels
if nargin<5
    l_data=max(length(positives),length(negatives));
    for label_ix=1:l_data
        graph_labels{label_ix}=['ROI',num2str(label_ix)];
    end
end



%% prepare cognitive domain based coloring
if color_mode==7
    % basic colouring is set here
    EmoInt=[255	192	0]/255;
    EmoInt2=[255	255	0]/255;
    MotVis=[0	32	96]/255;
    MotVis2=[0	112	192]/255;
    Visual=[0	255	204]/255;
    HiCog1=[112	48	160]/255;
    HiCog2=[131	60	12]/255;
    HiCog3=[0	102	0]/255;
    HiCog4=[255	80	80]/255;
    HiCog5=[190	0	0]/255;
    Artifa=[127 127 127]/255;
    Outsid=[0 0 0]/255;
    CogGrpCol=vertcat(EmoInt,MotVis,Visual,HiCog1,HiCog2,HiCog3,HiCog4,HiCog5,Artifa,Outsid);
    % at the moment colors are not connected to the atlases, but inferred
    % from the number of values submitted to the visualisation script
    if max(size(positives))==10 || max(size(positives))==11 || max(size(negatives))==10 || max(size(negatives))==11
        % SMITH10
        CogCol=[
            Visual
            Visual
            Visual
            HiCog1
            HiCog2
            MotVis
            HiCog4
            EmoInt
            HiCog3
            HiCog5
            Outsid
            ];
        
    elseif max(size(positives))==20 || max(size(positives))==21 || max(size(negatives))==20 || max(size(negatives))==21
        %BRAINMAP20
        CogCol=[
            EmoInt %  1
            EmoInt %  2
            EmoInt2 %  3 external
            EmoInt2 %  4 external
            EmoInt2 %  5 external
            MotVis2 %  6 visuospatial
            MotVis2 %  7 visuospatial
            MotVis %  8
            MotVis %  9
            Visual % 10
            Visual % 11
            Visual % 12
            HiCog1 % 13
            HiCog2 % 14
            HiCog3 % 15
            HiCog4 % 16
            HiCog4 % 17
            HiCog5 % 18
            Artifa % 19
            Artifa % 20
            Outsid % 21
            ];
    elseif max(size(positives))==70 || max(size(positives))==71 || max(size(negatives))==70 || max(size(negatives))==71
        %BRAINMAP70
        CogCol=[
            Visual %  1
            Visual %  2
            Visual %  3
            Visual %  4
            Visual %  5
            MotVis2 %  6
            MotVis2 %  7
            MotVis2 %  8
            MotVis2 %  9
            MotVis2 % 10
            HiCog3 % 11
            HiCog3 % 12
            HiCog3 % 13
            EmoInt % 14
            MotVis2 % 15
            MotVis2 % 16
            EmoInt % 17
            EmoInt % 18
            EmoInt % 19
            EmoInt % 20
            EmoInt % 21
            EmoInt % 22
            EmoInt % 23
            HiCog1 % 24
            EmoInt % 25
            EmoInt % 26
            HiCog3 % 27
            HiCog1 % 28
            HiCog1 % 29
            HiCog1 % 30
            EmoInt % 31
            Visual % 32
            MotVis % 33
            MotVis % 34
            MotVis % 35
            MotVis % 36
            MotVis % 37
            HiCog1 % 38
            EmoInt % 39
            MotVis % 40
            EmoInt % 41
            Visual % 42
            Visual % 43
            HiCog4 % 44
            HiCog4 % 45
            HiCog4 % 46
            HiCog4 % 47
            HiCog5 % 48
            HiCog5 % 49
            HiCog5 % 50
            HiCog5 % 51
            EmoInt2 % 52
            HiCog3 % 53
            EmoInt2 % 54
            EmoInt2 % 55
            EmoInt2 % 56
            EmoInt2 % 57
            EmoInt2 % 58
            HiCog2 % 59
            HiCog2 % 60
            Visual % 61
            HiCog1 % 62
            HiCog4 % 63
            Artifa % 64
            EmoInt2 % 65
            HiCog2 % 66
            Artifa % 67
            Artifa % 68
            Artifa % 69
            Artifa % 70
            Outsid % 71
            ];
    else % this  is a failsafe if the data is not SMITH10, BRAINMAP20 or BRAINMAP70-based
        color_mode=1;
    end
end


%% correction, this may cause some conflict and misinterpretations so it still needs to be checked
positives(find(positives<0))=NaN;

%% positive values for negatives enforced
if ~isempty(negatives)
    if length(find(negatives(1,:)<=0))>length(find(negatives(1,:)>0))
        negatives(1,:)=-negatives(1,:);
    end
end
negatives(find(negatives<0))=NaN;

%% arrange vectors
if size(positives,1)>size(positives,2)
    positives=positives';
end
if size(negatives,1)>size(negatives,2)
    negatives=negatives';
end

%% normalise to [0..1]
if ~isempty(positives)
    if size(positives,1)>1
        max_positives=max(positives(1,:)+positives(2,:));
    else
        max_positives=max(positives(1,:));
    end
else
    max_positives=0;
end
if ~isempty(negatives)
    if size(negatives,1)>1
        max_negatives=max(negatives(1,:)+negatives(2,:));
    else
        max_negatives=max(negatives(1,:));
    end
else
    max_negatives=0;
end
if ~forced_max
    max_data=max(max_positives,max_negatives);
    %forced_ticks=forced_ticks/max_data;
else
    max_data=forced_max;
end
%if max_data>1
positives=positives/max_data;
negatives=negatives/max_data;
%end

%% check for error or SD data
if size(positives,1)==2
    positives_err1=positives(1,:)+positives(2,:);
    positives_err2=positives(1,:)-positives(2,:);
    positives=positives(1,:);
    draw_pos_err=1;
else
    draw_pos_err=0;
end
if size(negatives,1)==2
    negatives_err1=negatives(1,:)+negatives(2,:);
    negatives_err2=negatives(1,:)-negatives(2,:);
    negatives=negatives(1,:);
    draw_neg_err=1;
else
    draw_neg_err=0;
end

%% set defaults
draw_rose=0;
draw_line=0;
draw_polar=0;
polar_filled=0;
polar_overlap=0;
mirror_type=0;
draw_bar=0;

switch lower(type_of_chart)
    case {'rose','rosec',1}
        draw_rose=1;
        type_of_chart=1;
    case {'roselr',2}
        mirror_type=1;
        draw_rose=1;
        type_of_chart=2;
    case {'line','linec',3}
        draw_line=1;
        type_of_chart=3;
    case {'linelr',4}
        mirror_type=1;
        draw_line=1;
        type_of_chart=4;
    case {'polar','polafc',5}
        draw_polar=1;
        type_of_chart=5;
    case {'polarlr',6}
        mirror_type=1;
        draw_polar=1;
        type_of_chart=6;
    case {'polaroverlap',7}
        draw_polar=1;
        polar_overlap=1;
        type_of_chart=7;
    case {'polarf','polarfc',8}
        draw_polar=1;
        polar_filled=1;
        type_of_chart=8;
    case {'polarflr',9}
        mirror_type=1;
        draw_polar=1;
        polar_filled=1;
        type_of_chart=9;
    case{'barv',10}
        draw_bar=1;
        type_of_chart=10;
    case{'barh',11}
        draw_bar=2;
        type_of_chart=11;
    case{'roseovr',12}
        polar_overlap=1;
        draw_rose=2;
        type_of_chart=12;
    otherwise
end

switch color_mode
    case {'rb','br',1}
        poscol=[1 0 0] ;
        negcol=[0 0 1] ;
    otherwise
        poscol=[0.75 0.75 .75];
        negcol=[0.25 0.25 .25];
end





h0=figure;


clf
hold on
if draw_bar
    l_data=max(length(positives),length(negatives));
    if draw_bar==1
        if ~isempty(positives)
            bar(positives*max_data,'facecolor',poscol,'barwidth',1,'linewidth',2)
            ylim_max=1.1*max_data;
        else
            ylim_max=0.1*max_data;
        end
        if ~isempty(negatives)
            bar(-negatives*max_data,'facecolor',negcol,'barwidth',1,'linewidth',2)
            ylim_min=-1.1*max_data;
        else
            ylim_min=-0.1*max_data;
        end
        ylim([ylim_min ylim_max])
        xlim([0 l_data+1])
        line([0 l_data+1],[ylim_max ylim_max],'color',[0.25, 0.25, 0.25],'linestyle','-')
        line([l_data+1 l_data+1],[ylim_min ylim_max],'color',[0.25, 0.25, 0.25],'linestyle','-')
        set(gca,'ytick',[-forced_ticks(end:-1:1)*max_data 0 forced_ticks*max_data])
        for t_ix=1:length(forced_ticks)
            plot([0 l_data+1],[forced_ticks(t_ix) forced_ticks(t_ix)]*max_data,'k:')
            plot([0 l_data+1],-[forced_ticks(t_ix) forced_ticks(t_ix)]*max_data,'k:')
        end
        if draw_pos_err
            for e_ix=1:length(positives_err1)
                line([e_ix e_ix],[positives_err1(e_ix)*max_data positives_err2(e_ix)*max_data],'LineWidth',4,'Color','k')
            end
        end
        if draw_neg_err
            for e_ix=1:length(negatives_err1)
                line([e_ix e_ix],[-negatives_err1(e_ix)*max_data -negatives_err2(e_ix)*max_data],'LineWidth',4,'Color','k')
            end
        end
        set(gca,'xtick',[1:l_data],'xticklabel',graph_labels)
        axis square
    else
        if ~isempty(positives)
            barh(positives*max_data,'facecolor',poscol,'barwidth',1,'linewidth',2)
            xlim_max=1.1*max_data;
        else
            xlim_max=0.1*max_data;
        end
        if ~isempty(negatives)
            barh(-negatives*max_data,'facecolor',negcol,'barwidth',1,'linewidth',2)
            xlim_min=-1.1*max_data;
        else
            xlim_min=-0.1*max_data;
        end
        xlim([xlim_min xlim_max])
        ylim([0 l_data+1])
        line([xlim_max xlim_max],[0 l_data+1],'Color',[0.25, 0.25, 0.25],'linestyle','-')
        line([xlim_min xlim_max],[0 0],'Color',[0.25, 0.25, 0.25],'linestyle','-')
        set(gca,'xtick',[-forced_ticks(end:-1:1)*max_data 0 forced_ticks*max_data])
        for t_ix=1:length(forced_ticks)
            plot([forced_ticks(t_ix) forced_ticks(t_ix)]*max_data,[0 l_data+1],'k:')
            plot(-[forced_ticks(t_ix) forced_ticks(t_ix)]*max_data,[0 l_data+1],'k:')
        end
        
        
        if draw_pos_err
            for e_ix=1:length(positives_err1)
                line([positives_err1(e_ix)*max_data positives_err2(e_ix)*max_data],[e_ix e_ix],'LineWidth',4,'Color','k')
            end
        end
        if draw_neg_err
            for e_ix=1:length(negatives_err1)
                line([-negatives_err1(e_ix)*max_data -negatives_err2(e_ix)*max_data],[e_ix e_ix],'LineWidth',4,'Color','k')
            end
        end
        set(gca,'ytick',[1:l_data],'yticklabel',graph_labels)
        axis square ij
    end
    if nargin>5
        title(graph_title,'Interpreter','none')
    end
else % polar representations
    hold on
    % outer rectangle
    rectangle('position',[-1.5 -1.5 3 3],'curvature',[0 0],'linestyle','-','facecolor','w','edgecolor','w')
    
    %% draw the coordinate system
    % outer circle
    rectangle('position',[-1 -1 2 2],'curvature',[1 1],'linestyle','-','facecolor','w','edgecolor','k','Linewidth',polarlinewidth)
    plot(0,0,'k+')
    
    %% draw the results
    if polar_overlap
        data_points=max(size(positives,2),size(negatives,2))+2;
    else
        data_points=(size(negatives,2)+size(positives,2))+2;
    end
    single_output=0;
    if isempty(negatives) | isempty(positives) | polar_overlap
        data_points=data_points-2;
        single_output=1;
    end
    data_angles=((2*pi):-(2*pi/data_points):(2*pi/data_points))+(pi/2)-(pi/data_points);
    data_angles=data_angles-(pi/data_points);
    %% triangle representation
    signs={'positives','negatives'};
    for sign_ix=1:length(signs)
        do_draw_polar=0;
        switch signs{sign_ix}
            case 'positives'
                if ~isempty(positives)
                    do_draw_polar=1;
                    if single_output
                        draw_angles=[data_angles data_angles(1)]+(2*pi/data_points);
                        line_angles=[data_angles data_angles(1)]+(pi/data_points);
                    else
                        draw_angles=[data_angles(1:(length(positives)+1))]+(pi/data_points);
                        line_angles=[data_angles data_angles(1)];
                    end
                    draw_radii=[positives positives(1)];
                    if draw_pos_err
                        draw_radii1=[positives_err1 positives_err1(1)];
                        draw_radii2=[positives_err2 positives_err2(1)];
                    else
                        draw_radii1=[positives positives(1)];
                        draw_radii2=[positives positives(1)];
                    end
                    
                end
            case 'negatives'
                if ~isempty(negatives)
                    do_draw_polar=1;
                    if single_output
                        draw_angles=[data_angles data_angles(1)]+(2*pi/data_points);
                        line_angles=[data_angles data_angles(1)]+(pi/data_points);
                    else
                        if mirror_type
                            draw_angles=[data_angles((end):-1:(end-length(negatives)-2))]+(pi/data_points);
                            line_angles=[data_angles((end):-1:(end-length(negatives)-2))]+(2*pi/data_points);
                        else
                            draw_angles=[data_angles((end-length(negatives)):1:(end))]+(pi/data_points);
                            line_angles=[data_angles((end-length(negatives)):1:(end))];%-(pi/data_points);
                        end
                    end
                    draw_radii=[negatives negatives(1)];
                    if draw_neg_err
                        draw_radii1=[negatives_err1 negatives_err1(1)];
                        draw_radii2=[negatives_err2 negatives_err2(1)];
                    else
                        draw_radii1=[negatives negatives(1)];
                        draw_radii2=[negatives negatives(1)];
                    end
                end
            otherwise
                do_draw_polar=0;
        end
        
        if do_draw_polar
            %% calculate
            clear patchx patchy linex liney polarx polary polarx1 polary1 polarx2 polary2 errlinex elliney
            patchx2=[];
            patchy2=[];
            for data_ix=1:(length(draw_radii)-1)
                [patchx{data_ix},patchy{data_ix}]=pol2cart(draw_angles([data_ix data_ix data_ix+1 data_ix+1]),[0 draw_radii([data_ix data_ix]) 0]);
                [patchx2tmp,patchy2tmp]=pol2cart(draw_angles([data_ix data_ix+1 data_ix+1]),[draw_radii([data_ix data_ix data_ix+1])]);
                patchx2=[patchx2,patchx2tmp];
                patchy2=[patchy2,patchy2tmp];
                [rosex{data_ix},rosey{data_ix}]=pol2cart(draw_angles([data_ix data_ix+1 data_ix+1]),[draw_radii([data_ix data_ix data_ix+1])]);
                [linex{data_ix},liney{data_ix}]=pol2cart(line_angles([data_ix data_ix]),[0 draw_radii(data_ix)]);
                [errlinex{data_ix},errliney{data_ix}]=pol2cart(line_angles([data_ix data_ix]),[draw_radii1(data_ix) draw_radii2(data_ix)]);
                [ticksx{data_ix},ticksy{data_ix}]=pol2cart(line_angles([data_ix data_ix]),[1 1.025]);
                [polarx(data_ix),polary(data_ix)]=pol2cart(line_angles(data_ix),draw_radii(data_ix));
                [polarx1(data_ix),polary1(data_ix)]=pol2cart(line_angles(data_ix),draw_radii1(data_ix));
                [polarx2(data_ix),polary2(data_ix)]=pol2cart(line_angles(data_ix),draw_radii2(data_ix));
                [textx{data_ix},texty{data_ix}]=pol2cart(line_angles(data_ix),1.125);
                switch type_of_chart
                    case {1,2}
                        [ticklinex{data_ix},tickliney{data_ix}]=pol2cart(draw_angles([data_ix data_ix]),[0 1]);
                        [ticklinex2{data_ix},tickliney2{data_ix}]=pol2cart(draw_angles([data_ix+1 data_ix+1]),[0 1]);
                    otherwise
                        [ticklinex{data_ix},tickliney{data_ix}]=pol2cart(line_angles([data_ix data_ix]),[0 1]);
                end
            end
            %% draw
            for data_ix=1:(length(draw_radii)-1)
                
                switch color_mode
                    case {1}
                        poscol='r';
                        negcol='b';
                        edgecol='k';
                    case {2}
                        poscol=[0.75 0.75 .75];
                        negcol=[0.25 0.25 .25];
                        edgecol='k';
                    case {3}
                        poscol=draw_radii(data_ix);
                        negcol=-draw_radii(data_ix);
                        colormap(jet)
                        if single_output
                            caxis([0 1])
                            locmap=jet(101);
                        else
                            caxis([-1 1])
                            locmap=jet(201);
                        end
                        edgecol='k';
                    case {4}
                        poscol=draw_radii(data_ix);
                        negcol=-draw_radii(data_ix);
                        colormap(parula)
                        locmap=jet(101);
                        if single_output
                            caxis([0 1])
                            locmap=parula(101);
                        else
                            caxis([-1 1])
                            locmap=parula(201);
                        end
                        edgecol='k';
                    case {5}
                        poscol=[1 0+draw_radii(data_ix)/1 0+draw_radii(data_ix)/3];
                        negcol=[0+draw_radii(data_ix)/2 0+draw_radii(data_ix)/1 1];
                        edgecol='k';
                    case {6}
                        poscol=[0.5+draw_radii(data_ix)/2 0+draw_radii(data_ix)/2 0+draw_radii(data_ix)/3];
                        negcol=[0+draw_radii(data_ix)/4 0+draw_radii(data_ix)/2 0.5+draw_radii(data_ix)/2];
                        edgecol='k';
                    case {7}
                        poscol=CogCol(data_ix,:);
                        negcol=CogCol(data_ix,:);
                        edgecol=[0.5 0.5 0.5];
                    otherwise
                        poscol='r';
                        negcol='b';
                        edgecol='k';
                end
                switch signs{sign_ix}
                    case 'positives'
                        drawcol=poscol;
                    case 'negatives'
                        drawcol=negcol;
                end
                if draw_rose
                    if draw_rose==2 && strcmp(signs{sign_ix},'negatives')
                        line(rosex{data_ix},rosey{data_ix},'Color',[0.25, 0.25, 0.25],'LineWidth',2)
                    else
                        patch(patchx{data_ix},patchy{data_ix},drawcol,'facealpha',0.5,'Edgecolor',edgecol)
                    end
                    if draw_pos_err || draw_neg_err
                        line(errlinex{data_ix},errliney{data_ix},'LineWidth',2,'Color',[0.25, 0.25, 0.25],'LineStyle','-')
                    end
                end
                if draw_line
                    switch color_mode
                        case{3,4}
                            if single_output
                                poscol=locmap(round(draw_radii(data_ix)*100)+1,:);
                                negcol=locmap(round((1-draw_radii(data_ix))*100)+1,:);
                            else
                                poscol=locmap(round(draw_radii(data_ix)*100)+100,:);
                                negcol=locmap(round((1-draw_radii(data_ix))*100)+1,:);
                            end
                        otherwise
                    end
                    switch signs{sign_ix}
                        case 'positives'
                            drawcol=poscol;
                        case 'negatives'
                            drawcol=negcol;
                    end
                    line(linex{data_ix},liney{data_ix},'LineWidth',4,'Color',drawcol)
                end
                if draw_polar
                    switch color_mode
                        case{3,4}
                            poscol=locmap(end,:);
                            negcol=locmap(1,:);
                            switch signs{sign_ix}
                                case 'positives'
                                    drawcol=poscol;
                                case 'negatives'
                                    drawcol=negcol;
                            end
                        otherwise
                    end
                    if single_output
                        polarx=[polarx polarx(1)];
                        polary=[polary polary(1)];
                        polarx1=[polarx1 polarx1(1)];
                        polary1=[polary1 polary1(1)];
                        polarx2=[polarx2 polarx2(1)];
                        polary2=[polary2 polary2(1)];
                    else
                        polarx=[polarx 0 polarx(1)];
                        polary=[polary 0 polary(1)];
                        polarx1=[polarx1 0 polarx1(1)];
                        polary1=[polary1 0 polary1(1)];
                        polarx2=[polarx2 0 polarx2(1)];
                        polary2=[polary2 0 polary2(1)];
                    end
                    if polar_filled
                        patch(polarx,polary,drawcol,'facealpha',0.5)
                    else
                        patch(polarx,polary,drawcol,'facecolor','none','edgecolor',drawcol,'linewidth',1,'edgealpha',0.2)
                        if draw_pos_err || draw_neg_err
                            patch(polarx1,polary1,drawcol,'facecolor','none','edgecolor',drawcol,'linewidth',1,'linestyle',':','edgealpha',0.2)
                            patch(polarx2,polary2,drawcol,'facecolor','none','edgecolor',drawcol,'linewidth',1,'linestyle',':','edgealpha',0.2)
                        end
                    end
                    
                    line(linex{data_ix},liney{data_ix},'LineWidth',1,'Color',[0.25, 0.25, 0.25],'linestyle','-')
                end
                

                %% ticks
                line(ticksx{data_ix},ticksy{data_ix},'Linewidth',polarlinewidth*2,'Color',[0.25, 0.25, 0.25],'linestyle','-')
                line(ticklinex{data_ix},tickliney{data_ix},'Linewidth',polarlinewidth,'Color',[0.25, 0.25, 0.25],'linestyle',':')
                if type_of_chart==1 | type_of_chart==2 & data_ix==(length(draw_radii)-1)
                    line(ticklinex2{data_ix},tickliney2{data_ix},'Linewidth',polarlinewidth,'Color',[0.25, 0.25, 0.25],'linestyle',':')
                end
                if polar_text_perpendicular
                    if mod((line_angles(data_ix)/(2*pi)*360)+90,360)>180
                        text(textx{data_ix},texty{data_ix},graph_labels{data_ix},'horizontalalignment','right','verticalalignment','middle','rotation',180+line_angles(data_ix)/(2*pi)*360,'interpreter','none')
                    else
                        text(textx{data_ix},texty{data_ix},graph_labels{data_ix},'horizontalalignment','left','verticalalignment','middle','rotation',line_angles(data_ix)/(2*pi)*360,'interpreter','none')
                    end
                else
                    if ~isempty(findstr(graph_labels{data_ix},'BM'))
                        text(textx{data_ix},texty{data_ix},regexprep(graph_labels{data_ix},'[^\w'']',''),'horizontalalignment','center','verticalalignment','bottom','rotation',-90+line_angles(data_ix)/(2*pi)*360, 'interpreter','none', 'FontSize',9)
                    else
                        text(textx{data_ix},texty{data_ix},graph_labels{data_ix},'horizontalalignment','center','verticalalignment','bottom','rotation',-90+line_angles(data_ix)/(2*pi)*360,'interpreter','none','FontSize',9)
                    end
                end
            end

            
            if nargin>5
                title(graph_title,'Interpreter','none')
            end
            if ~single_output 
                if mirror_type
                    line([0 0],[-1.025 1.025],'Linewidth',polarlinewidth*2,'Color',[0.25, 0.25, 0.25],'linestyle','-')
                else
                    line([0 0],[-1.025 1.025],'Linewidth',polarlinewidth,'Color',[0.25, 0.25, 0.25],'linestyle','-')
                end
            else
                line([0 0],[-1.025 1.025],'Linewidth',polarlinewidth,'Color',[0.25, 0.25, 0.25],'linestyle',':')
            end
            
            %% draw the ticklines and ticklabels
            ticks_at=forced_ticks;
            for t_ix=1:length(ticks_at)
                rectangle('position',[-ticks_at(t_ix) -ticks_at(t_ix) 2*ticks_at(t_ix) 2*ticks_at(t_ix)],'curvature',[1 1],'linestyle',':','edgecolor',[0.25, 0.25, 0.25],'Linewidth',polarlinewidth)
                text('Interpreter','none','String',[' ',num2str(ticks_at(t_ix)*max_data)],'Position',[0 ticks_at(t_ix)],'horizontalalignment','left','verticalalignment','bottom','FontSize',8)
                text('Interpreter','none','String',[' ',num2str(ticks_at(t_ix)*max_data)],'Position',[0 -ticks_at(t_ix)],'horizontalalignment','left','verticalalignment','top','interpreter','none','FontSize',8)
            end
            axis([-1.5 2 -1.5 2])
            axis equal
            axis off
        end
    end
end
