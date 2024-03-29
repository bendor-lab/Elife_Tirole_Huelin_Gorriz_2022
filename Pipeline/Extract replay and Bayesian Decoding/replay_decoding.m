% REPLAY DECODING
% Extracts spikes for each replay event and decodes it. Saves in a
% structure called repay_track. Inside each field (each track) the replay events are saved.
% Loads: extracted_replay_events, extracted_clusters and extracted_place_fields_BAYESIAN

function decoded_replay_events = replay_decoding(varargin)

    % Load parameters
    parameters = list_of_parameters;
    load('extracted_replay_events.mat');
    load('extracted_clusters.mat');
    load('extracted_place_fields_BAYESIAN.mat');

    if ~isempty(varargin)
        switch varargin{1}
            case 'control_fixed_rate'
                option= 2;
            case 'control_fixed_count'
                option= 3;
            case 'replay_rate_shuffle'
                 option= 4;
            otherwise
                option= 1;
        end
    else
        option= 1;
    end
    
    % REPLAY EVENTS STRUCTURE
    %replay_events is an empty template for replay event analysis.
    %Each track will create its own field

    replay_events = struct('replay_id',{},...%the id of the candidate replay events in chronological order
        'spikes',{});%column 1 is spike id, column 2 is spike time

    % TAKE SPIKES FROM ONLY good place fields (on at least one track)
    
    sorted_spikes = zeros(size(clusters.spike_id));
    sorted_spikes(:,1) = clusters.spike_id;
    sorted_spikes(:,2) = clusters.spike_times;
    
    all_units = unique(clusters.spike_id);
    non_pyramidal = setdiff(all_units,place_fields_BAYESIAN.good_place_cells);
    for i = 1 : length(non_pyramidal)
        non_pyramidal_indices = find(sorted_spikes(:,1)== non_pyramidal(i));
        sorted_spikes(non_pyramidal_indices,:) = [];
        non_pyramidal_indices =[];
    end
    num_spikes = length(sorted_spikes);
    num_units = length(place_fields_BAYESIAN.good_place_cells);


    % EXTRACT SPIKES IN REPLAY EVENTS

    num_replay = size(replay.onset, 2);
    current_replay = 1;
    current_replay_spikes = [];
    
    for i = 1 : num_spikes
        if sorted_spikes(i,2) > replay.offset(current_replay)
            replay_events(current_replay).replay_id = current_replay;
            replay_events(current_replay).spikes = current_replay_spikes;
            current_replay = current_replay + 1;
            if current_replay > num_replay
                break
            end
            current_replay_spikes = [];
        end
        if sorted_spikes(i,2) >= replay.onset(current_replay)
            current_replay_spikes = [current_replay_spikes; sorted_spikes(i,:)];
        end
    end
 
    num_replay_events = length(replay_events);
    msg = [num2str(num_replay_events), ' candidate events.'];
    disp(msg);

    % Save all replay events all tracks
    for j = 1:length(place_fields_BAYESIAN.track)
        decoded_replay_events(j).replay_events = replay_events;
    end

    %%%%%% BAYESIAN DECODING ON REPLAY EVENTS %%%%%%
    replay_starts = replay.onset;
    replay_ends = replay.offset;
    
    % Get time vectors for bayesian decoding and matrix with spike count
    if option== 1 % normal replay decoding
        disp('Spike count...');
        replayEvents_bayesian_spike_count = spike_count([],replay_starts,replay_ends,'N');  % Takes a long time to run (few minutes)
        save replayEvents_bayesian_spike_count replayEvents_bayesian_spike_count
        
        % Run bayesian decoding
        disp('Decoding position...');
        estimated_position = bayesian_decoding([],'replayEvents_bayesian_spike_count','N');
    elseif option== 2 % fixed rate control
        disp('Fixed spike rate control...');
        modify_spike_rate_count;
        disp('Decoding position control...');
        estimated_position = bayesian_decoding([],'replayEvents_FIXED_spike','N');
    elseif option== 3 % fixed count control
        disp('Fixed spike count control...');
        modify_spike_count;
        disp('Decoding position control...');
        estimated_position = bayesian_decoding([],'replayEvents_FIXED_spike','N');
    elseif option == 4 % replay rate shuffle
        disp('replay rate shuffle');
        shuffle_spike_rate_replay;
        estimated_position = bayesian_decoding([],'replayEvents_bayesian_spike_count','N');
    end
  
    % Save in structure
    for j = 1:length(place_fields_BAYESIAN.track)
        for i = 1 : num_replay_events
            decoded_replay_events(j).replay_events(i).timebins_edges = estimated_position(j).replay_events(i).replay_time_edges; 
            decoded_replay_events(j).replay_events(i).timebins_centre = estimated_position(j).replay_events(i).replay_time_centered; 
            decoded_replay_events(j).replay_events(i).timebins_index = 1:length(estimated_position(j).replay_events(i).replay_time_centered);
            decoded_replay_events(j).replay_events(i).decoded_position = estimated_position(j).replay_events(i).replay; % normalized by all tracks   
        end
    end
      
    % Saves structure
    if option== 1 || option == 4% normal replay decoding
        save decoded_replay_events decoded_replay_events
    elseif option== 2 || option ==3 % fixed rate/count controls
        save decoded_replay_events_FIXED decoded_replay_events
    end
    
end


