% SPEARMAN WITHOUT BAYESIAN BIAS
% Spearman correlation on replay events without using bayesian bias to select between multi track events
% (events significant for more than one track). Creates new data control folders
% INPUT - spearman_modification:
% no_multievents: re-run analysis pipeline having removed the multi track events from the spearman structure
% pvalue_selection: instead of selecting multi track events with bayesian bias, uses the smaller p-value to assign them to a track.
% Then re-runs analysis pipeline

function spearman_without_BayesianBias(spearman_modification)

% Load data structures
if exist('significant_replay_events_spearman_individual_exposures.mat','file') % for re-rexposure sessions
    load significant_replay_events_spearman_individual_exposures.mat
    rexp = 1;
else
    load significant_replay_events_spearman.mat
    rexp = [];
end
load extracted_position
load extracted_sleep_state
load extracted_place_fields_BAYESIAN

current_folder = strsplit(pwd,'\');

% Remove multi events for tracks
if strcmp(spearman_modification,'no_multievents')
    
    multi_events = significant_replay_events.multi_tracks_index;
    
    % Find multi events indices in each track and deletes the relevant events
    for t = 1 : length(significant_replay_events.track)
     
        [~,multi_events_idxs,~] = intersect(significant_replay_events.track(t).index,multi_events); %finds indices for multi events
        
        significant_replay_events.track(t).index(multi_events_idxs) = [];
        significant_replay_events.track(t).ref_index(multi_events_idxs) = [];
        significant_replay_events.track(t).event_times(multi_events_idxs) = [];
        significant_replay_events.track(t).replay_score(multi_events_idxs) = [];
        significant_replay_events.track(t).p_value(multi_events_idxs) = [];
        significant_replay_events.track(t).bayesian_bias(multi_events_idxs) = [];
        significant_replay_events.track(t).event_segment_best_score(multi_events_idxs) = [];
        significant_replay_events.track(t).spikes(multi_events_idxs) = [];
        significant_replay_events.track(t).decoded_position(multi_events_idxs) = [];
        significant_replay_events.track(t).event_duration(multi_events_idxs) = [];
    end
    
    cd('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\CONTROLS\Spearman_no_multiEvents')
    
    % Create a new folder with the name of the data analysed
    if exist(current_folder{end})~=7
        mkdir(current_folder{end})
    end
    % Go to control data folder
    cd(['X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\CONTROLS\Spearman_no_multiEvents\' current_folder{end}])
    % Save structure needed for analysis pipeline
    if ~isempty(rexp)
        save significant_replay_events_spearman_individual_exposures.mat significant_replay_events
    else
        save significant_replay_events_spearman.mat significant_replay_events
    end
    save extracted_position.mat position
    save extracted_sleep_state.mat sleep_state
    save extracted_place_fields_BAYESIAN.mat place_fields_BAYESIAN
    
    
    % Run again analysis
    if isempty(rexp)
        % no reexposure
        sort_replay_events([],'spearman');
    else
        % reexposure
        sort_replay_events([1,3;2,4],'spearman');
    end
    
    % individual folders
    rate_remapping_TRACK_PAIRS([],'spearman');
    plot_rate_remapping('TRACK_PAIRS','spearman');
    
    save_path = 'X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\CONTROLS\Figures\Spearman_no_multiEvents';
    savefig([save_path '\' current_folder{end} '_spearman.fig'])
    saveas(gcf,[save_path '\' current_folder{end} '_spearman.png'])
    close all
    
% Take multi events for tracks and classify based on p-value   
elseif strcmp(spearman_modification,'pvalue_selection')
    
    load scored_replay.mat
    load scored_replay_segments.mat
    multi_events = significant_replay_events.multi_tracks_index;
    
    % Find multi events indices in each track and deletes the relevant events
    for t = 1 : 2 %just look at the first two tracks for now
     
        [~,multi_events_idxs,~] = intersect(significant_replay_events.track(t).index,multi_events); %finds indices for multi events
        ref_index = significant_replay_events.track(t).ref_index(multi_events_idxs);% trace back reference index
        
        % For each multi-event find the p-value that the event had in the other track
        for e = 1 : length(ref_index)
            T1_pvals = min([scored_replay(1).replay_events(ref_index(e)).spearman_p scored_replay1(1).replay_events(ref_index(e)).spearman_p scored_replay2(1).replay_events(ref_index(e)).spearman_p]); 
            T2_pvals = min([scored_replay(2).replay_events(ref_index(e)).spearman_p scored_replay1(2).replay_events(ref_index(e)).spearman_p scored_replay2(2).replay_events(ref_index(e)).spearman_p]); 
            [~,sig_track] = min([T1_pvals T2_pvals]);
            % If the new significant track is not the same than the current
            % one, remove it from the list and change it to the pertinent track
            if sig_track ~= t 
                if t == 1
                    other_t = 2;
                else
                    other_t = 1;
                end
                significant_replay_events.track(other_t).index(end+1) = significant_replay_events.track(t).index(multi_events_idxs(e)); % copy info to the other track
                significant_replay_events.track(t).index(multi_events_idxs(e)) = []; % delete from current track
                significant_replay_events.track(other_t).ref_index(end+1) = significant_replay_events.track(t).ref_index(multi_events_idxs(e));
                significant_replay_events.track(t).ref_index(multi_events_idxs(e)) = [];
                significant_replay_events.track(other_t).event_times(end+1) = significant_replay_events.track(t).event_times(multi_events_idxs(e));
                significant_replay_events.track(t).event_times(multi_events_idxs(e)) = [];
                significant_replay_events.track(other_t).replay_score(end+1) = significant_replay_events.track(t).replay_score(multi_events_idxs(e));
                significant_replay_events.track(t).replay_score(multi_events_idxs(e)) = [];
                significant_replay_events.track(other_t).p_value(end+1) = significant_replay_events.track(t).p_value(multi_events_idxs(e));
                significant_replay_events.track(t).p_value(multi_events_idxs(e)) = [];
                significant_replay_events.track(other_t).bayesian_bias(end+1) = significant_replay_events.track(t).bayesian_bias(multi_events_idxs(e));
                significant_replay_events.track(t).bayesian_bias(multi_events_idxs(e)) = [];
                significant_replay_events.track(other_t).event_segment_best_score(end+1) = significant_replay_events.track(t).event_segment_best_score(multi_events_idxs(e));
                significant_replay_events.track(t).event_segment_best_score(multi_events_idxs(e)) = [];
                significant_replay_events.track(other_t).spikes(end+1) = significant_replay_events.track(t).spikes(multi_events_idxs(e));
                significant_replay_events.track(t).spikes(multi_events_idxs(e)) = [];
                significant_replay_events.track(other_t).decoded_position(end+1) = significant_replay_events.track(t).decoded_position(multi_events_idxs(e));
                significant_replay_events.track(t).decoded_position(multi_events_idxs(e)) = [];
                significant_replay_events.track(other_t).event_duration(end+1) = significant_replay_events.track(t).event_duration(multi_events_idxs(e));
                significant_replay_events.track(t).event_duration(multi_events_idxs(e)) = [];
            end
    end
    
    cd('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\CONTROLS\Spearman_pvalue_selection')
    
    % Create a new folder with the name of the data analysed
    if exist(current_folder{end})~=7
        mkdir(current_folder{end})
    end
    cd(['X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\CONTROLS\Spearman_pvalue_selection\' current_folder{end}])
     % Save structure needed for analysis pipeline
    if ~isempty(rexp)
        save significant_replay_events_spearman_individual_exposures.mat significant_replay_events
    else
        save significant_replay_events_spearman.mat significant_replay_events
    end
    save extracted_position.mat position
    save extracted_sleep_state.mat sleep_state
    save extracted_place_fields_BAYESIAN.mat place_fields_BAYESIAN
    
    % Run again analysis
      if isempty(rexp)
        % no reexposure
        sort_replay_events([],'spearman');
    else
        % reexposure
        sort_replay_events([1,3;2,4],'spearman');
    end
    
    % individual folders
    rate_remapping_TRACK_PAIRS([],'spearman');
    plot_rate_remapping('TRACK_PAIRS','spearman');
    
    save_path = 'X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\CONTROLS\Figures\Spearman_pvalue_selection';
    savefig([save_path '\' current_folder{end} '_spearman.fig'])
    saveas(gcf,[save_path '\' current_folder{end} '_spearman.png'])
    close all
     
end


end