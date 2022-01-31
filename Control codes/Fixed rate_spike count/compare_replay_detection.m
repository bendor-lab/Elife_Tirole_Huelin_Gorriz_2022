function compare_replay_detection()
% load significant events
load('sorted_replay_wcorr_FIXED.mat')
% load all candidate events
sorted_replay_fixed= sorted_replay;
load('scored_replay_FIXED.mat');
load('significant_replay_events_wcorr_FIXED.mat');

curr_folder= pwd;
delim= strfind(curr_folder, '\');
session_name= curr_folder(delim(end)+1:end);

% load original sorted replay
load(['..\..\..\' session_name '\sorted_replay_wcorr.mat']);

% now look at how they were classified
for i=1:length(significant_replay_events.pre_ripple_threshold_index) % all candidate events
    if ismember(i,sorted_replay(1).index.sleepPOST) && ismember(i,sorted_replay_fixed(1).index.sleepPOST) % detected as track 1 for both methods
        event_coord(i,:)= [i i];
        event_cdata(i,:)= [1 0 0]; % red for track 1
        event_marker{i}= '.';
    elseif ismember(i,sorted_replay(2).index.sleepPOST) && ismember(i,sorted_replay_fixed(2).index.sleepPOST) % detected as track 2 for both methods
        event_coord(i,:)= [i i];
        event_cdata(i,:)= [0 0 1]; % blue for track 2
        event_marker{i}= '.';
    elseif ismember(i,sorted_replay(1).index.sleepPOST) && ismember(i,sorted_replay_fixed(2).index.sleepPOST) % detected as opposite tracks by the two methods
        event_coord(i,:)= [i-10 i+10]; % offset to control
        event_cdata(i,:)= [1 0 0]; % should be track 1
        event_marker{i}= '*';  % asterix for misclassified by control
    elseif ismember(i,sorted_replay(2).index.sleepPOST) && ismember(i,sorted_replay_fixed(1).index.sleepPOST) % detected as opposite tracks by the two methods
        event_coord(i,:)= [i+10 i-10]; % offset to normal
        event_cdata(i,:)= [0 0 1]; % should be track 2
        event_marker{i}= '*';  % asterix for misclassified by control
    elseif (ismember(i,sorted_replay(1).index.sleepPOST) && ~ismember(i,sorted_replay_fixed(1).index.sleepPOST))  % event only detected by normal method
        event_coord(i,:)= [i+10 i-10]; % offset to normal
        event_cdata(i,:)= [1 0 0]; % should be track 1
        event_marker{i}= 'x';  % cross for missing in control
    elseif (ismember(i,sorted_replay(2).index.sleepPOST) && ~ismember(i,sorted_replay_fixed(2).index.sleepPOST))
        event_coord(i,:)= [i+10 i-10]; % offset to normal
        event_cdata(i,:)= [0 0 1]; % should be track 2
        event_marker{i}= 'x';  % cross for missing in control
    elseif (~ismember(i,sorted_replay(1).index.sleepPOST) && ismember(i,sorted_replay_fixed(1).index.sleepPOST))  % event only detected by control method
        event_coord(i,:)= [i-10 i+10]; % offset to control
        event_cdata(i,:)= [1 0 0]; % should be track 1
        event_marker{i}= 'd';  % diamond for missing in normal
    elseif (~ismember(i,sorted_replay(2).index.sleepPOST) && ismember(i,sorted_replay_fixed(2).index.sleepPOST))
        event_coord(i,:)= [i-10 i+10]; % offset to control
        event_cdata(i,:)= [0 0 1]; % should be track 2
        event_marker{i}= 'd';  % diamond for missing in normal
    elseif (~ismember(i,sorted_replay(1).index.sleepPOST) && ~ismember(i,sorted_replay_fixed(1).index.sleepPOST)) &&... % not a significant event
           (~ismember(i,sorted_replay(2).index.sleepPOST) && ~ismember(i,sorted_replay_fixed(2).index.sleepPOST))
        event_coord(i,:)= [i i]; % offset to control
        event_cdata(i,:)= [0.5 0.5 0.5]; % neither tracks
        event_marker{i}= '.';  % 
    end
end

min_index= min([sorted_replay(1).index.sleepPOST sorted_replay(2).index.sleepPOST sorted_replay_fixed(1).index.sleepPOST sorted_replay_fixed(2).index.sleepPOST])
figure
hold on;
arrayfun(@(i) plot(event_coord(i,1),event_coord(i,2),'Marker',event_marker{i},'Color',event_cdata(i,:)),min_index:length(event_coord));
xlabel('standard')
ylabel('control')
title({'Red - track1 , Blue- track2, dots - similar both methods,'; 'diamond - detected in control only, Cross- detected in standard only'})
end