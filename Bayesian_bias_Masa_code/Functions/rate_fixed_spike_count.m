function rate_fixed_spike_count(bayesian_spike_count,timebin,option)

% load real spike count
if strcmp(bayesian_spike_count,'replayEvents_common_good_cell_bayesian_spike_count')
    load('replayEvents_common_good_cell_bayesian_spike_count.mat')
    replayEvents_bayesian_spike_count = replayEvents_common_good_cell_bayesian_spike_count;
else
    load('replayEvents_bayesian_spike_count.mat');
end

parameters=list_of_parameters;

% those variables are unchanged
replayEvents_FIXED_spike_rate.replay_events_indices= replayEvents_bayesian_spike_count.replay_events_indices;
replayEvents_FIXED_spike_rate.replay_events= replayEvents_bayesian_spike_count.replay_events;

replayEvents_FIXED_spike_count.replay_events_indices = replayEvents_bayesian_spike_count.replay_events_indices;
replayEvents_FIXED_spike_count.replay_events = replayEvents_bayesian_spike_count.replay_events;

% allocate empty array
replayEvents_FIXED_spike_rate.n.replay= NaN(size(replayEvents_bayesian_spike_count.n.replay));

num_events= length(unique(replayEvents_bayesian_spike_count.replay_events_indices(~isnan(replayEvents_bayesian_spike_count.replay_events_indices))));


cd 'rate_fixed'

if timebin < 1 % if smaller than 1, then it is specifying ms instead of number of time bins
        % find average FR cell during events (regardless of significance)
        cell_spikes_events= zeros(size(replayEvents_bayesian_spike_count.n.replay,1),1);
        for this_event=1:num_events
            event_idx= (replayEvents_bayesian_spike_count.replay_events_indices == this_event);
            
            % real spikes count
            spike_count_real= replayEvents_bayesian_spike_count.n.replay(:,event_idx);
            
            % keep only cells where there was at least one spike
            num_spikes_per_cell= sum(spike_count_real,2);
            non_zero_spikes= find(num_spikes_per_cell>0);
            
            % before do the mean look if it has a value above 0
            zero_idx= cell_spikes_events(non_zero_spikes)==0;
            
            % count number of spikes for first timers
            cell_spikes_events(non_zero_spikes(zero_idx))= num_spikes_per_cell(non_zero_spikes(zero_idx));
            % calculate average number of spikes for the rest
            cell_spikes_events(non_zero_spikes(~zero_idx))= mean([cell_spikes_events(non_zero_spikes(~zero_idx)) num_spikes_per_cell(non_zero_spikes(~zero_idx))],2);
        end
        
        % it feels stupid doing it again but...
        for this_event=1:num_events
            event_idx= (replayEvents_bayesian_spike_count.replay_events_indices == this_event);
            spike_count_real= replayEvents_bayesian_spike_count.n.replay(:,event_idx);
            num_spikes_per_cell= sum(spike_count_real,2);
            non_zero_spikes= find(num_spikes_per_cell>0);
            spike_num_fixed= zeros(size(spike_count_real));
            
            % scale number of spikes by average FR
            spikes_this_event= num_spikes_per_cell(non_zero_spikes);
            spike_num_fixed(non_zero_spikes,:)= replayEvents_bayesian_spike_count.n.replay(non_zero_spikes,event_idx).*(cell_spikes_events(non_zero_spikes)./spikes_this_event);
            
            replayEvents_FIXED_spike_count.n.replay(:,event_idx)= spike_num_fixed;
        end
        
        save('replayEvents_FIXED_spike_count_20ms.mat','replayEvents_FIXED_spike_count');
    
elseif timebin == 1
    if strcmp(option,'median spike fix')
        % find Median FR cell during events (regardless of significance)
        cell_spikes_events= zeros(size(replayEvents_bayesian_spike_count.n.replay,1),1);
        for ncell = 1:length(zeros(size(replayEvents_bayesian_spike_count.n.replay,1),1))
            k = 1;
            for this_event=1:num_events
                event_idx= (replayEvents_bayesian_spike_count.replay_events_indices == this_event);
                % real spikes count
                spike_count_real= replayEvents_bayesian_spike_count.n.replay(:,event_idx);
                
                if spike_count_real(ncell) ~= 0 % if the cell fires during this event
                    spike_count_this_cell(k) = spike_count_real(ncell);
                    k = k + 1;
                end
            end
            cell_spikes_events(ncell) = median(spike_count_this_cell);% calculate the median spike count for this cell
        end
        
        % it feels stupid doing it again but...
        for this_event=1:num_events
            event_idx= (replayEvents_bayesian_spike_count.replay_events_indices == this_event);
            spike_count_real= replayEvents_bayesian_spike_count.n.replay(:,event_idx);
            num_spikes_per_cell= sum(spike_count_real,2);
            non_zero_spikes= find(num_spikes_per_cell>0);
            spike_num_fixed= zeros(size(spike_count_real));
            
            % fix at median spike count
            spike_num_fixed(non_zero_spikes,:)= cell_spikes_events(non_zero_spikes);
            
            replayEvents_FIXED_spike_count.n.replay(:,event_idx)= spike_num_fixed;
        end
                
        save('replayEvents_FIXED_spike_count_one_bin_median_spike_count.mat','replayEvents_FIXED_spike_count');
      
    elseif isempty(option)
        
        % find average FR cell during events (regardless of significance)
        cell_spikes_events= zeros(size(replayEvents_bayesian_spike_count.n.replay,1),1);
        for this_event=1:num_events
            event_idx= (replayEvents_bayesian_spike_count.replay_events_indices == this_event);
            
            % real spikes count
            spike_count_real= replayEvents_bayesian_spike_count.n.replay(:,event_idx);
            
            % keep only cells where there was at least one spike
            num_spikes_per_cell= sum(spike_count_real,2);
            non_zero_spikes= find(num_spikes_per_cell>0);
            
            % before do the mean look if it has a value above 0
            zero_idx= cell_spikes_events(non_zero_spikes)==0;
            
            % count number of spikes for first timers
            cell_spikes_events(non_zero_spikes(zero_idx))= num_spikes_per_cell(non_zero_spikes(zero_idx));
            % calculate average number of spikes for the rest
            cell_spikes_events(non_zero_spikes(~zero_idx))= mean([cell_spikes_events(non_zero_spikes(~zero_idx)) num_spikes_per_cell(non_zero_spikes(~zero_idx))],2);
        end
        
        % it feels stupid doing it again but...
        for this_event=1:num_events
            event_idx= (replayEvents_bayesian_spike_count.replay_events_indices == this_event);
            spike_count_real= replayEvents_bayesian_spike_count.n.replay(:,event_idx);
            num_spikes_per_cell= sum(spike_count_real,2);
            non_zero_spikes= find(num_spikes_per_cell>0);
            spike_num_fixed= zeros(size(spike_count_real));
            
            % scale number of spikes by average FR
            spikes_this_event= num_spikes_per_cell(non_zero_spikes);
            spike_num_fixed(non_zero_spikes,:)= replayEvents_bayesian_spike_count.n.replay(non_zero_spikes,event_idx).*(cell_spikes_events(non_zero_spikes)./spikes_this_event);
            
            replayEvents_FIXED_spike_count.n.replay(:,event_idx)= spike_num_fixed;
        end
        
        save('replayEvents_FIXED_spike_count_one_bin.mat','replayEvents_FIXED_spike_count');
    end
end

cd ..
end