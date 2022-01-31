function modify_spike_rate_count()
% load real spike count
load('replayEvents_bayesian_spike_count.mat');
parameters=list_of_parameters;

% those variables are unchanged
replayEvents_FIXED_spike_rate.replay_events_indices= replayEvents_bayesian_spike_count.replay_events_indices;
replayEvents_FIXED_spike_rate.replay_events= replayEvents_bayesian_spike_count.replay_events;
% allocate empty array
replayEvents_FIXED_spike_rate.n.replay= NaN(size(replayEvents_bayesian_spike_count.n.replay));

num_events= length(unique(replayEvents_bayesian_spike_count.replay_events_indices(~isnan(replayEvents_bayesian_spike_count.replay_events_indices))));

% find average FR cell during events (regardless of significance)
cell_FR_events= zeros(size(replayEvents_bayesian_spike_count.n.replay,1),1);
for this_event=1:num_events
    event_idx= find(replayEvents_bayesian_spike_count.replay_events_indices == this_event);
    
    % real spikes count
    spike_count_real= replayEvents_bayesian_spike_count.n.replay(:,event_idx);
    
    % keep only events where there was at least one spike
    num_spikes_per_cell= sum(spike_count_real,2);
    non_zero_spikes= find(num_spikes_per_cell>0);
    zero_idx= cell_FR_events(non_zero_spikes)==0;
    
    % count number of spikes
    cell_FR_events(non_zero_spikes(zero_idx))= num_spikes_per_cell(non_zero_spikes(zero_idx))/(length(event_idx)*parameters.replay_bin_width);
    % calculate average firing rate during event
    cell_FR_events(non_zero_spikes(~zero_idx))= mean([cell_FR_events(non_zero_spikes(~zero_idx)) num_spikes_per_cell(non_zero_spikes(~zero_idx))/(length(event_idx)*parameters.replay_bin_width)],2);   
end

% it feels stupid doing it again but...
for this_event=1:num_events
    event_idx= find(replayEvents_bayesian_spike_count.replay_events_indices == this_event);
    spike_count_real= replayEvents_bayesian_spike_count.n.replay(:,event_idx);  
    num_spikes_per_cell= sum(spike_count_real,2);
    non_zero_spikes= find(num_spikes_per_cell>0);
    spike_FR_fixed= zeros(size(spike_count_real));
    
    % scale number of spikes by average FR
    FR_this_event= num_spikes_per_cell(non_zero_spikes)/(length(event_idx)*parameters.replay_bin_width);
    spike_FR_fixed(non_zero_spikes,:)= replayEvents_bayesian_spike_count.n.replay(non_zero_spikes,event_idx).*(cell_FR_events(non_zero_spikes)./FR_this_event);
    
    replayEvents_FIXED_spike_rate.n.replay(:,event_idx)= spike_FR_fixed;
end

save('replayEvents_FIXED_spike_rate.mat','replayEvents_FIXED_spike_rate');

end