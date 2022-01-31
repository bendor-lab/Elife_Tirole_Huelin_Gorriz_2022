function control_detection_FIXED_spike_rate(current_folder,RExp)
% works in separate folder (CONTROLS\fixed_rate)
% input list of folders to process
% and vector of reexposure options

if exist('..\CONTROLS\fixed_rate')~=7
    mkdir('..\CONTROLS\fixed_rate');
end

    
    if exist(['..\CONTROLS\fixed_rate\' current_folder])~=7
        mkdir(['..\CONTROLS\fixed_rate\' current_folder]);
    end
    
    % copy over files needed 
    copyfile('extracted_replay_events.mat',['..\CONTROLS\fixed_rate\' current_folder]);
    copyfile('decoded_replay_events.mat',['..\CONTROLS\fixed_rate\' current_folder]);
    copyfile('decoded_replay_events_segments.mat',['..\CONTROLS\fixed_rate\' current_folder]);
    copyfile('extracted_clusters.mat',['..\CONTROLS\fixed_rate\' current_folder]);
    copyfile('extracted_place_fields_BAYESIAN.mat',['..\CONTROLS\fixed_rate\' current_folder]);
    copyfile('extracted_position.mat',['..\CONTROLS\fixed_rate\' current_folder]);
    copyfile('replayEvents_bayesian_spike_count.mat',['..\CONTROLS\fixed_rate\' current_folder]);
    copyfile('extracted_sleep_state.mat',['..\CONTROLS\fixed_rate\' current_folder]);

    cd(['..\CONTROLS\fixed_rate\' current_folder]);    
    
    replay_decoding('control_fixed_rate');
    
    disp('scoring replay events')
    load('decoded_replay_events_FIXED.mat');
    scored_replay = replay_scoring(decoded_replay_events,[0 1 0 0]); % weighted corr  
    save scored_replay_FIXED scored_replay;

    % RUN SHUFFLES
    tic;
    disp('running shuffles')
    num_shuffles=1000;
    analysis_type=[0 1 0 0];  %just weighted correlation 
    p = gcp; % Starting new parallel pool
    shuffle_choice={'PRE spike_train_circular_shift','PRE place_field_circular_shift', 'POST place bin circular shift'};
    tic
    if ~isempty(p)
        for shuffle_id=1:length(shuffle_choice)
            shuffle_type{shuffle_id}.shuffled_track = parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events,'control_fixed_spike');
        end
    else
        disp('parallel processing not possible');
        for shuffle_id=1:length(shuffle_choice)
            shuffle_type{shuffle_id}.shuffled_track = run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events,'control_fixed_spike');
        end
    end
    save shuffled_tracks_FIXED shuffle_type;
    disp('time to run shuffles was...');
    toc
    % Evaluate significance
    load('scored_replay_FIXED.mat');
    load('shuffled_tracks.mat');
    scored_replay= replay_significance(scored_replay, shuffle_type);
    save scored_replay_FIXED scored_replay

    %%%%%%analyze segments%%%%%%%%%%
    % splitting replay events
    tic
    replay_decoding_split_events('control_fixed_spike');
    load decoded_replay_events_segments_FIXED;
    scored_replay1 = replay_scoring(decoded_replay_events1,[0 1 0 0]);
    scored_replay2 = replay_scoring(decoded_replay_events2,[0 1 0 0]);
    save scored_replay_segments_FIXED scored_replay1 scored_replay2;
    num_shuffles=1000;
    analysis_type=[0 1 0 0];  %just weighted correlation

    load decoded_replay_events_segments_FIXED;
    p = gcp; % Starting new parallel pool
    if ~isempty(p)
        for shuffle_id=1:length(shuffle_choice)
            shuffle_type1{shuffle_id}.shuffled_track = parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events1,'control_fixed_spike');
            shuffle_type2{shuffle_id}.shuffled_track = parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events2,'control_fixed_spike');
        end
    else
        disp('parallel processing not possible');
        for shuffle_id=1:length(shuffle_choice)
            shuffle_type1{shuffle_id}.shuffled_track = run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events1,'control_fixed_spike');
            shuffle_type2{shuffle_id}.shuffled_track = run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events2,'control_fixed_spike');
        end
    end

    save shuffled_tracks_segments_FIXED shuffle_type1 shuffle_type2;

    load scored_replay_segments_FIXED; load shuffled_tracks_segments_FIXED;
    scored_replay1=replay_significance(scored_replay1, shuffle_type1);
    scored_replay2=replay_significance(scored_replay2, shuffle_type2);
    save scored_replay_segments_FIXED scored_replay1 scored_replay2
    disp('running segments took...')
    toc

    %%%% analyze significant replays from whole events and segmented events
    % only using replay events passing threshold for ripple power
  
    if isempty(RExp) % no reexposure
        number_of_significant_replays(0.05,3,'wcorr',[],'control_fixed_spike'); % pval, ripple zscore, method, reexposure
    else
        % reexposure
        number_of_significant_replays(0.05,3,'wcorr',2,'control_fixed_spike');
    end
    
    if isempty(RExp) % no reexposure
        sort_replay_events([],'control_fixed_spike'); 
    else
        % reexposure
        sort_replay_events([1 3;2 4],'control_fixed_spike'); 
    end
    
    rate_remapping_TRACK_PAIRS([],'control_fixed_spike');
    plot_rate_remapping('TRACK_PAIRS','control_fixed_rate');
    
%     compare_replay_detection;

end