% Compare replay events detected by different methods
% INPUT:
    % period: sleep, awake, or [] for both awake+sleep
    % varargin: 'control_fixed_spike' or 'rate_detection_control'
    
function compare_replay_detection(period,varargin)

if isempty(varargin)
    error('Missing control input')
elseif ~isempty(varargin) && strcmp(varargin{1},'control_fixed_spike')
    % load significant events
    load('sorted_replay_wcorr_FIXED.mat')
    % load all candidate events
    sorted_replay_control = sorted_replay;
    load('scored_replay_FIXED.mat');
    load('significant_replay_events_wcorr_FIXED.mat');
elseif ~isempty(varargin) && strcmp(varargin{1},'rate_detection_control')
    % load significant events
    load('sorted_replay_wcorr_RATE.mat')
    % load all candidate events
    sorted_replay_control = sorted_replay;
    load('scored_replay.mat');
    if exist('significant_replay_events_wcorr_individual_exposures.mat','file')
        load('significant_replay_events_wcorr_individual_exposures.mat');
    else
        load('significant_replay_events_wcorr.mat');
    end
end

curr_folder = pwd;
delim= strfind(curr_folder, '\');
session_name= curr_folder(delim(end)+1:end);

% load original sorted replay
load(['..\..\..\' session_name '\sorted_replay_wcorr.mat']);
t1_count = 0;
t2_count = 0;
all_t1 = 0;
all_t2 = 0;

% now look at how they were classified
for i = 1 : length(significant_replay_events.pre_ripple_threshold_index) % all candidate events
    if ismember(i,sorted_replay(1).index.(strcat(period,'POST'))) && ismember(i,sorted_replay_control(1).index.(strcat(period,'POST'))) % detected as track 1 for both methods
        event_coord(i,:)= [i i];
        event_cdata(i,:)= [1 0 0]; % red for track 1
        event_marker{i}= '.';
        t1_count = t1_count +1;
        all_t1 = all_t1 +1;
    elseif ismember(i,sorted_replay(2).index.(strcat(period,'POST'))) && ismember(i,sorted_replay_control(2).index.(strcat(period,'POST'))) % detected as track 2 for both methods
        event_coord(i,:)= [i i];
        event_cdata(i,:)= [0 0 1]; % blue for track 2
        event_marker{i}= '.';
        t2_count = t2_count +1;
        all_t2 = all_t2 +1;
    elseif ismember(i,sorted_replay(1).index.(strcat(period,'POST'))) && ismember(i,sorted_replay_control(2).index.(strcat(period,'POST'))) % detected as opposite tracks by the two methods
        event_coord(i,:)= [i-10 i+10]; % offset to control
        event_cdata(i,:)= [1 0 0]; % should be track 1
        event_marker{i}= '*';  % red asterix for misclassified by control
        all_t1 = all_t1 +1;
    elseif ismember(i,sorted_replay(2).index.(strcat(period,'POST'))) && ismember(i,sorted_replay_control(1).index.(strcat(period,'POST'))) % detected as opposite tracks by the two methods
        event_coord(i,:)= [i+10 i-10]; % offset to normal
        event_cdata(i,:)= [0 0 1]; % should be track 2
        event_marker{i}= '*';  % blue asterix for misclassified by control
        all_t2 = all_t2 +1;
    elseif (ismember(i,sorted_replay(1).index.(strcat(period,'POST'))) && ~ismember(i,sorted_replay_control(1).index.(strcat(period,'POST'))))  % event only detected by normal method
        event_coord(i,:)= [i+10 i-10]; % offset to normal
        event_cdata(i,:)= [1 0 0]; % should be track 1
        event_marker{i}= 'x';  % red cross for missing in control
        all_t1 = all_t1 +1;
    elseif (ismember(i,sorted_replay(2).index.(strcat(period,'POST'))) && ~ismember(i,sorted_replay_control(2).index.(strcat(period,'POST'))))
        event_coord(i,:)= [i+10 i-10]; % offset to normal
        event_cdata(i,:)= [0 0 1]; % should be track 2
        event_marker{i}= 'x';  % blue cross for missing in control
        all_t2 = all_t2 +1;
    elseif (~ismember(i,sorted_replay(1).index.(strcat(period,'POST'))) && ismember(i,sorted_replay_control(1).index.(strcat(period,'POST'))))  % event only detected by control method
        event_coord(i,:)= [i-10 i+10]; % offset to control
        event_cdata(i,:)= [1 0 0]; % should be track 1
        event_marker{i}= 'd';  % red diamond for missing in normal
    elseif (~ismember(i,sorted_replay(2).index.(strcat(period,'POST'))) && ismember(i,sorted_replay_control(2).index.(strcat(period,'POST'))))
        event_coord(i,:)= [i-10 i+10]; % offset to control
        event_cdata(i,:)= [0 0 1]; % should be track 2
        event_marker{i}= 'd';  % blue diamond for missing in normal
    elseif (~ismember(i,sorted_replay(1).index.(strcat(period,'POST'))) && ~ismember(i,sorted_replay_control(1).index.(strcat(period,'POST')))) &&... % not a significant event
           (~ismember(i,sorted_replay(2).index.(strcat(period,'POST'))) && ~ismember(i,sorted_replay_control(2).index.(strcat(period,'POST'))))
        event_coord(i,:)= [i i]; % offset to control
        event_cdata(i,:)= [0.5 0.5 0.5]; % grey neither tracks
        event_marker{i}= '.';  % 
    end
end

min_index= min([sorted_replay(1).index.(strcat(period,'POST')) sorted_replay(2).index.(strcat(period,'POST')) sorted_replay_control(1).index.(strcat(period,'POST')) sorted_replay_control(2).index.(strcat(period,'POST'))]);

%Calculate proportion of overlap detection:
total_num_sig_events_POST = length(sorted_replay(1).index.POST)+length(sorted_replay(2).index.POST);
total_num_sig_events_POST_control = length(sorted_replay_control(1).index.POST)+length(sorted_replay_control(2).index.POST);
perc_overlap_num_sig_events_POST = total_num_sig_events_POST_control/total_num_sig_events_POST*100;
perc_overlap_T1_events = t1_count/all_t1*100;
perc_overlap_T2_events = t2_count/all_t2*100;

figure
hold on;
arrayfun(@(i) plot(event_coord(i,1),event_coord(i,2),'Marker',event_marker{i},'Color',event_cdata(i,:)),min_index:length(event_coord));
xlabel('standard')
ylabel('control')
annotation('textbox',[.2 .8 .1 .1],'String',{'Red - track1 , Blue- track2, dots - similar both methods,'; 'diamond - detected in control only, Cross- detected in standard only'})
title({['% overlap num sig POST events = ' num2str(perc_overlap_num_sig_events_POST)],[' % overlap num sig T1 events = ' num2str(perc_overlap_T1_events)] ...
    [' % overlap num sig T2 events = ' num2str(perc_overlap_T2_events)]})
end