function shuffle_spike_rate_replay(varargin)

if isempty(varargin)
    method= 'between_tracks';
elseif strcmp(varargin{1},'intrinsic_bias')
    method= 'within_tracks';
end

% load real spike count
load('replayEvents_bayesian_spike_count.mat');
load('extracted_place_fields_BAYESIAN.mat');
load('decoded_replay_events.mat');
if exist('significant_replay_events_wcorr.mat') ==2
    load('significant_replay_events_wcorr.mat');
elseif exist('significant_replay_events_wcorr_individual_exposures.mat') ==2
    load('significant_replay_events_wcorr_individual_exposures.mat');
end
parameters=list_of_parameters;

% define kernel
n_sigma= 0.1; % 50ms
gauss_ker= gausswin(n_sigma/(1/1000));
gauss_ker= gauss_ker./sum(gauss_ker); % one isolated spike should have an amplitude of 1

% those variables are unchanged
replayEvents_SCALED_rate.replay_events_indices= replayEvents_bayesian_spike_count.replay_events_indices;
replayEvents_SCALED_rate.replay_events= replayEvents_bayesian_spike_count.replay_events;
% allocate empty array
replayEvents_SCALED_rate.n.replay= NaN(size(replayEvents_bayesian_spike_count.n.replay));

for this_track=1:2
    % find significant events
    sig_event_track_idx= significant_replay_events.track(this_track).ref_index;
    % find average FR cell during events (regardless of significance)
    cell_FR_events{this_track}= zeros(size(replayEvents_bayesian_spike_count.n.replay,1),1); % length number of cells
for this_event=1:length(sig_event_track_idx)
    event_idx= find(replayEvents_bayesian_spike_count.replay_events_indices == sig_event_track_idx(this_event));
    
    % real spikes count
    replay_time= significant_replay_events.track(this_track).event_times(this_event);
    replay_dur= significant_replay_events.track(this_track).event_duration(this_event);
    replay_spikes= significant_replay_events.track(this_track).spikes{this_event};
    ts_event_edges= min(replay_time-(replay_dur/2),min(replay_spikes(:,2)))-n_sigma :1/1000: max(replay_time+(replay_dur/2),max(replay_spikes(:,2)))+n_sigma;
    
    % find cells that were active 
    [active_cells,cell_event_idx,cell_idx]= intersect(significant_replay_events.track(this_track).spikes{this_event}(:,1),place_fields_BAYESIAN.good_place_cells);
    for this_cell=1:length(active_cells)
        binned_spike_train= histcounts(replay_spikes(replay_spikes(:,1) == cell_event_idx(this_cell),2),ts_event_edges);
        if cell_FR_events{this_track}(cell_idx(this_cell)) ~=0
            cell_FR_events{this_track}(cell_idx(this_cell))= mean([cell_FR_events{this_track}(cell_idx(this_cell)) max(filter(gauss_ker,1,binned_spike_train))]);
        else
            cell_FR_events{this_track}(cell_idx(this_cell))= max(filter(gauss_ker,1,binned_spike_train));
        end    
    end

end
end

num_events= length(unique(replayEvents_bayesian_spike_count.replay_events_indices(~isnan(replayEvents_bayesian_spike_count.replay_events_indices))));
 
% assign scaled rate to each event regardless of significance, from
% distribution 
switch method
    case 'between_tracks'
        rate_distribution= sort(vertcat(cell_FR_events{:}));
        rate_distribution(rate_distribution==0)=[];
    case 'within_tracks'
        
end
fraction_of_cells = (0:(length(rate_distribution)-1))/(length(rate_distribution)-1);
num_events= length(unique(replayEvents_bayesian_spike_count.replay_events_indices(~isnan(replayEvents_bayesian_spike_count.replay_events_indices))));
for this_event=1:num_events
    event_idx= find(replayEvents_bayesian_spike_count.replay_events_indices == this_event);
    random_rate= []; replay_rate_shuffled=[]; scaling_factor= [];
    
    replay_time= decoded_replay_events(1).replay_events(this_event).timebins_edges;
    replay_dur= replay_time(end)- replay_time(1);
    replay_spikes=decoded_replay_events(1).replay_events(this_event).spikes;
    ts_event_edges= min(replay_time-(replay_dur/2),min(replay_spikes(:,2)))-n_sigma :1/1000: max(replay_time+(replay_dur/2),max(replay_spikes(:,2)))+n_sigma;

    for this_cell=1:size(replayEvents_bayesian_spike_count.n.replay,1)
        binned_spike_train= histcounts(replay_spikes(replay_spikes(:,1) == place_fields_BAYESIAN.good_place_cells(this_cell),2),ts_event_edges);
        max_inst_FR= max(filter(gauss_ker,1,binned_spike_train));

        random_rate(this_cell)=  interp1(fraction_of_cells,rate_distribution,rand(1),'linear');

        % scale number of spikes by max inst FR 
        if max_inst_FR>0
            replay_rate_shuffled(this_cell,:)= replayEvents_bayesian_spike_count.n.replay(this_cell,event_idx).*(random_rate(this_cell)./max_inst_FR);
            scaling_factor(this_cell)= random_rate(this_cell)./max_inst_FR;
        else
            replay_rate_shuffled(this_cell,:)= replayEvents_bayesian_spike_count.n.replay(this_cell,event_idx);
            scaling_factor(this_cell)= 1; % no scaling
        end
    end
    
    replayEvents_SCALED_rate.n.replay(:,event_idx)= replay_rate_shuffled;
    replayEvents_SCALED_rate.scaling_factors(:,this_event)= scaling_factor;
    
end

replayEvents_bayesian_spike_count= replayEvents_SCALED_rate;
save('replayEvents_bayesian_spike_count.mat','replayEvents_bayesian_spike_count');

end