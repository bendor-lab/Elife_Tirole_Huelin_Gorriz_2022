function remapping = all_cells_per_track(spikes_control,folders,varargin)
% folders is [] for individual session, or cell array of folders
% varargin can be 'spearman' or 'wcorr'
if isempty(folders)
    a=pwd;
    folders={a};
    cd ..
end
cd(folders{1});

if ~isempty(varargin)
    method= varargin{1};
    disp(['running for... ' varargin{1}])
else
    method= 'wcorr';
end

% Compare properties of common cells between tracks during replay
for curr_track = 1 : 2
    data(curr_track).epoch = compare_replay_across_tracks(method,curr_track);

% Find common cells that are active from PRE to POST and save
[data(curr_track).epoch(1).PRE_to_POST_active_cells, data(curr_track).epoch(2).PRE_to_POST_active_cells ]= deal(intersect(data(curr_track).epoch(1).ID_active_cells_during_replay,data(curr_track).epoch(2).ID_active_cells_during_replay));
end 

% Allocate variables
for curr_track=1:size(data,2)
    for epoch=1:size(data(curr_track).epoch,2)
        remapping(curr_track).epoch(epoch).experiment=[];
        remapping(curr_track).epoch(epoch).raw_peak_BAYESIAN_plfield_1=[];
        remapping(curr_track).epoch(epoch).track1_mean_replay_spikes=[];
        remapping(curr_track).epoch(epoch).track1_mean_replay_spikes_nonZero=[];
        remapping(curr_track).epoch(epoch).track1_mean_replay_rate=[];
        
        remapping(curr_track).epoch(epoch).ID_active_cells_during_replay=[];
        remapping(curr_track).epoch(epoch).fraction_of_active_cells_during_replay=[];
        remapping(curr_track).epoch(epoch).PRE_to_POST_active_cells=[];
    end
end
%cd ..


for i=1:length(folders)
    cd(folders{i});
    %data = compare_replay_across_tracks(method,curr_track);
    % Find common cells that are active from PRE to POST and save
    %[data(1).PRE_to_POST_active_cells, data(2).PRE_to_POST_active_cells ]= deal(intersect(data(2).ID_active_cells_during_replay,data(1).ID_active_cells_during_replay));
    
    for curr_track=1:size(data,2)
        for epoch=1:size(data(curr_track).epoch,2)
            remapping(curr_track).epoch(epoch).folder=folders{i};
            remapping(curr_track).epoch(epoch).experiment=[remapping(curr_track).epoch(epoch).experiment; (i*ones(size(data(curr_track).epoch(epoch).raw_peak_BAYESIAN_plfield_1)))];
            
            remapping(curr_track).epoch(epoch).raw_peak_BAYESIAN_plfield_1=[remapping(curr_track).epoch(epoch).raw_peak_BAYESIAN_plfield_1 data(curr_track).epoch(epoch).raw_peak_BAYESIAN_plfield_1];
            remapping(curr_track).epoch(epoch).track1_mean_replay_spikes=[remapping(curr_track).epoch(epoch).track1_mean_replay_spikes; data(curr_track).epoch(epoch).track1_mean_replay_spikes];
            remapping(curr_track).epoch(epoch).track1_mean_replay_rate=[remapping(curr_track).epoch(epoch).track1_mean_replay_rate; data(curr_track).epoch(epoch).track1_mean_replay_rate];
            remapping(curr_track).epoch(epoch).track1_mean_replay_spikes_nonZero=[remapping(curr_track).epoch(epoch).track1_mean_replay_spikes_nonZero; data(curr_track).epoch(epoch).track1_mean_replay_spikes_nonZero];
            
            remapping(curr_track).epoch(epoch).ID_active_cells_during_replay=[remapping(curr_track).epoch(epoch).ID_active_cells_during_replay; data(curr_track).epoch(epoch).ID_active_cells_during_replay];
            remapping(curr_track).epoch(epoch).fraction_of_active_cells_during_replay=[remapping(curr_track).epoch(epoch).fraction_of_active_cells_during_replay; data(curr_track).epoch(epoch).fraction_of_active_cells_during_replay];
            remapping(curr_track).epoch(epoch).PRE_to_POST_active_cells=[remapping(curr_track).epoch(epoch).PRE_to_POST_active_cells; data(curr_track).epoch(epoch).PRE_to_POST_active_cells];
            
        end
    end
    % cd ..
end

curr_folder = pwd;
fname = strsplit(a,'\');
cd(['X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\' fname{end}])

switch method
    case 'wcorr'
        save all_cells_analysis_single_track_wcorr_control remapping
    case 'spearman'
        save all_cells_analysis_single_track_spearman_control remapping
    otherwise
        save all_cells_analysis_single_track remapping
end
cd(curr_folder)

end


function remapping= compare_replay_across_tracks(method,curr_track)

switch method
    case 'wcorr'
        if exist('significant_replay_events_wcorr.mat')==2
            load('significant_replay_events_wcorr.mat');
        elseif exist('significant_replay_events_wcorr_individual_exposures.mat')==2
            load('significant_replay_events_wcorr_individual_exposures');
        end
        load('sorted_replay_wcorr');
        %significant_replay_events= significant_events_wcorr;
    case 'spearman'
        if exist('significant_replay_events_spearman.mat')==2
            load('significant_replay_events_spearman.mat');
        elseif exist('significant_replay_events_spearman_individual_exposures.mat')==2
            load('significant_replay_events_spearman_individual_exposures.mat');
        end
        load('sorted_replay_spearman');
%         significant_replay_events= significant_events_spearman;
    otherwise
        load('significant_replay_events');
        load('sorted_replay');
end

load('extracted_place_fields_BAYESIAN');

number_of_tracks=length(sorted_replay);

for epoch = 1:2  %PRE or POST
    if epoch==1
        index1  = sorted_replay(curr_track).index.sleepPRE; %indices replay events PRE for T1
        events1 = significant_replay_events.track(curr_track); %info from events in T1
    elseif epoch==2
        index1  = sorted_replay(curr_track).index.sleepPOST; %indices replay events POST for T1
        events1 = significant_replay_events.track(curr_track);
    end
    
    remapping(epoch).good_cells = [];
    remapping(epoch).fraction_of_active_cells_during_replay =[];
    remapping(epoch).track1_mean_replay_spikes=[];
    remapping(epoch).track1_mean_replay_spikes_nonZero=[];
    remapping(epoch).track1_mean_replay_rate=[];
    remapping(epoch).raw_peak_BAYESIAN_plfield_1=[];
    
    % Find BAYESIAN good cells in track
    good_cells = [place_fields_BAYESIAN.track(curr_track).good_cells];
    if spikes_control == 1
        peakFR = find(place_fields_BAYESIAN.track(curr_track).raw_peak(good_cells) >= 2);
        %2 spks/min
        remapping(epoch).good_cells = [place_fields_BAYESIAN.track(curr_track).good_cells(peakFR)];
    else
        remapping(epoch).good_cells = [place_fields_BAYESIAN.track(curr_track).good_cells];
    end
        
    
    % Calculate number of spikes in replay and divide by duration (based on 20ms binning, not exact time due to binning of spikes)
    if length(index1)>0  %if there's events for this track
        for j=1:length(index1) %for each event in T1
            for i=1:length(remapping(epoch).good_cells) %for each good cell
                % Find # spikes from good cell in this event
                remapping(epoch).spikes1(i,j) = length(find(events1.spikes{index1(j)}(:,1)==remapping(epoch).good_cells(i)));
                % Find FR from good cell in this event
                remapping(epoch).rate1(i,j) = length(find(events1.spikes{index1(j)}(:,1)==remapping(epoch).good_cells(i)))/events1.event_duration(index1(j));
            end
        end
        
        %%%% Track 1 - track 2 comparison of only good place cells (on both tracks)
        
        % Get place fields
        remapping(epoch).raw_peak_BAYESIAN_plfield_1 = place_fields_BAYESIAN.track(curr_track).raw_peak(remapping(epoch).good_cells);
        % Mean number of spikes across all events for each cell
        remapping(epoch).track1_mean_replay_spikes = mean(remapping(epoch).spikes1,2);
        % Sum of spikes per cell across all events divided by number of events where cell is active
        remapping(epoch).track1_mean_replay_spikes_nonZero = sum(remapping(epoch).spikes1,2)./sum(sign(remapping(epoch).spikes1),2); %mean of events with 1 or more spikes
        % Mean FR per cell across all events
        remapping(epoch).track1_mean_replay_rate = mean(remapping(epoch).rate1,2);
        
        % Find cells that have a mean FR >0 during replay & save cell IDs
        index = find(remapping(epoch).track1_mean_replay_spikes~=0);  %only analyze neurons that are active during replay.
        remapping(epoch).ID_active_cells_during_replay = remapping(epoch).good_cells(index);
        
        % Fraction of cells from good cells that are active in replay
        remapping(epoch).fraction_of_active_cells_during_replay = length(index)/length(remapping(epoch).track1_mean_replay_spikes);
        % Remove silent cells
        remapping(epoch).track1_mean_replay_rate = remapping(epoch).track1_mean_replay_rate(index);
        remapping(epoch).track1_mean_replay_spikes = remapping(epoch).track1_mean_replay_spikes(index);
        remapping(epoch).raw_peak_BAYESIAN_plfield_1 = remapping(epoch).raw_peak_BAYESIAN_plfield_1(index);
        remapping(epoch).track1_mean_replay_spikes_nonZero = remapping(epoch).track1_mean_replay_spikes_nonZero(index);
        
    end
end



% switch method
%     case 'wcorr'
%         save all_cells_analysis_TRACK_PAIRS_wcorr remapping
%     case 'spearman'
%         save all_cells_analysis_TRACK_PAIRS_spearman remapping
%     otherwise
%         save all_cells_analysis_TRACK_PAIRS remapping
% end

end