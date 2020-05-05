function remapping= rate_remapping_TRACK_PAIRS(folders,varargin)
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

data=compare_replay_across_tracks(method);
for epoch=1:size(data,1)
    for track_pair=1:size(data,2)
        remapping(epoch,track_pair).experiment=[];
        remapping(epoch,track_pair).place_field_diff=[];
        remapping(epoch,track_pair).replay_spike_diff=[];
        remapping(epoch,track_pair).replay_spike_diff_nonZero=[];
        remapping(epoch,track_pair).replay_rate_diff=[];
        remapping(epoch,track_pair).place_field_centre_diff=[];
        remapping(epoch,track_pair).place_fields_BAYESIAN1=[];
        remapping(epoch,track_pair).place_fields_BAYESIAN2=[];
        remapping(epoch,track_pair).track1_mean_replay_spikes=[];
        remapping(epoch,track_pair).track2_mean_replay_spikes=[];
        remapping(epoch,track_pair).track1_mean_replay_spikes_nonZero=[];
        remapping(epoch,track_pair).track2_mean_replay_spikes_nonZero=[];
        remapping(epoch,track_pair).track1_mean_replay_rate=[];
        remapping(epoch,track_pair).track2_mean_replay_rate=[];
    end
end
cd ..


for i=1:length(folders)
    cd(folders{i});
    data=compare_replay_across_tracks(method);
    
    for epoch=1:size(data,1)
        for track_pair=1:size(data,2)
            remapping(epoch,track_pair).folder=folders{i};
            remapping(epoch,track_pair).experiment=[remapping(epoch,track_pair).experiment; (i*ones(size(data(epoch,track_pair).place_field_diff)))];
            remapping(epoch,track_pair).place_field_diff=[remapping(epoch,track_pair).place_field_diff; data(epoch,track_pair).place_field_diff];
            remapping(epoch,track_pair).replay_spike_diff=[remapping(epoch,track_pair).replay_spike_diff; data(epoch,track_pair).replay_spike_diff];
             remapping(epoch,track_pair).replay_spike_diff_nonZero=[remapping(epoch,track_pair).replay_spike_diff_nonZero; data(epoch,track_pair).replay_spike_diff_nonZero];
         
            remapping(epoch,track_pair).replay_rate_diff=[remapping(epoch,track_pair).replay_rate_diff; data(epoch,track_pair).replay_rate_diff];
            remapping(epoch,track_pair).place_field_centre_diff=[remapping(epoch,track_pair).place_field_centre_diff; data(epoch,track_pair).place_field_centre_diff];
            remapping(epoch,track_pair).place_fields_BAYESIAN1=[remapping(epoch,track_pair).place_fields_BAYESIAN1 data(epoch,track_pair).place_fields_BAYESIAN1];
            remapping(epoch,track_pair).place_fields_BAYESIAN2=[remapping(epoch,track_pair).place_fields_BAYESIAN2 data(epoch,track_pair).place_fields_BAYESIAN2];
            remapping(epoch,track_pair).track1_mean_replay_spikes=[remapping(epoch,track_pair).track1_mean_replay_spikes; data(epoch,track_pair).track1_mean_replay_spikes];
            remapping(epoch,track_pair).track2_mean_replay_spikes=[remapping(epoch,track_pair).track2_mean_replay_spikes; data(epoch,track_pair).track2_mean_replay_spikes];
            remapping(epoch,track_pair).track1_mean_replay_rate=[remapping(epoch,track_pair).track1_mean_replay_rate; data(epoch,track_pair).track1_mean_replay_rate];
            remapping(epoch,track_pair).track2_mean_replay_rate=[remapping(epoch,track_pair).track2_mean_replay_rate; data(epoch,track_pair).track2_mean_replay_rate];
        end
    end
    cd ..
end

switch method
    case 'wcorr'
        save rate_remapping_analysis_TRACK_PAIRS_wcorr remapping
    case 'spearman'
        save rate_remapping_analysis_TRACK_PAIRS_spearman remapping
    otherwise
        save rate_remapping_analysis_TRACK_PAIRS remapping
end

end


function remapping= compare_replay_across_tracks(method)

switch method
    case 'wcorr'
        load('significant_replay_events_wcorr');
        load('sorted_replay_wcorr');
%         significant_replay_events= significant_events_wcorr;
    case 'spearman'
        load('significant_replay_events_spearman');
        load('sorted_replay_spearman');
%         significant_replay_events= significant_events_spearman;
    otherwise
        load('significant_replay_events');
        load('sorted_replay');
end

load('extracted_place_fields_BAYESIAN');

number_of_tracks=length(sorted_replay);

if number_of_tracks==2
    track_combinations=[1 2];
elseif number_of_tracks==3
    track_combinations=[1 2; 2 3; 3 1];
elseif number_of_tracks==4
    track_combinations=[1 2; 2 3; 3 4; 4 1; 1 3; 2 4];
end

for track_pair=1:size(track_combinations,1)
    tracks_to_compare=track_combinations(track_pair,:);
    for epoch=1:2  %PRE or POST
        if epoch==1
            index1=sorted_replay(tracks_to_compare(1)).index.sleepPRE;
            events1=significant_replay_events.track(tracks_to_compare(1));
            index2=sorted_replay(tracks_to_compare(2)).index.sleepPRE;
            events2=significant_replay_events.track(tracks_to_compare(2));
        elseif epoch==2
            index1=sorted_replay(tracks_to_compare(1)).index.sleepPOST;
            events1=significant_replay_events.track(tracks_to_compare(1));
            index2=sorted_replay(tracks_to_compare(2)).index.sleepPOST;
            events2=significant_replay_events.track(tracks_to_compare(2));
        end
        remapping(epoch,track_pair).fraction_of_cells_not_silent_during_replay=[];
        remapping(epoch,track_pair).place_field_diff=[];
        remapping(epoch,track_pair).replay_spike_diff=[];
        remapping(epoch,track_pair).replay_rate_diff=[];
        remapping(epoch,track_pair).place_field_centre_diff=[];
        remapping(epoch,track_pair).tracks_compared=[];
        remapping(epoch,track_pair).track1_mean_replay_spikes=[];
        remapping(epoch,track_pair).track2_mean_replay_spikes=[];
        remapping(epoch,track_pair).track1_mean_replay_rate=[];
        remapping(epoch,track_pair).track2_mean_replay_rate=[];
        remapping(epoch,track_pair).place_fields_BAYESIAN1=[];
        remapping(epoch,track_pair).place_fields_BAYESIAN2=[];     
        remapping(epoch,track_pair).good_cells=intersect(place_fields_BAYESIAN.track(tracks_to_compare(1)).good_cells, place_fields_BAYESIAN.track(tracks_to_compare(2)).good_cells);
        
        if length(index1)>0 & length(index2)>0
            for j=1:length(index1)
                for i=1:length(remapping(epoch,track_pair).good_cells)
                    remapping(epoch,track_pair).spikes1(i,j)=length(find(events1.spikes{index1(j)}(:,1)==remapping(epoch,track_pair).good_cells(i)));
                    remapping(epoch,track_pair).rate1(i,j)=length(find(events1.spikes{index1(j)}(:,1)==remapping(epoch,track_pair).good_cells(i)))/events1.event_duration(index1(j));
                end
            end
            for j=1:length(index2)
                for i=1:length(remapping(epoch,track_pair).good_cells) %calculate number of spikes in replay and divide by duration (based on binning, not exact time due to binning of spikes)
                    remapping(epoch,track_pair).spikes2(i,j)=length(find(events2.spikes{index2(j)}(:,1)==remapping(epoch,track_pair).good_cells(i)));
                    remapping(epoch,track_pair).rate2(i,j)=length(find(events2.spikes{index2(j)}(:,1)==remapping(epoch,track_pair).good_cells(i)))/events2.event_duration(index2(j));
                end
            end
            
            %track 1 - track 2 comparison of only good place cells (on both tracks)
            remapping(epoch,track_pair).place_fields_BAYESIAN1=place_fields_BAYESIAN.track(tracks_to_compare(1)).raw_peak(remapping(epoch,track_pair).good_cells);
            remapping(epoch,track_pair).place_fields_BAYESIAN2=place_fields_BAYESIAN.track(tracks_to_compare(2)).raw_peak(remapping(epoch,track_pair).good_cells);
            
            remapping(epoch,track_pair).place_field_centre_diff=(abs(place_fields_BAYESIAN.track(tracks_to_compare(1)).centre_of_mass(remapping(epoch,track_pair).good_cells)-...
                place_fields_BAYESIAN.track(tracks_to_compare(2)).centre_of_mass(remapping(epoch,track_pair).good_cells)))';
            
            remapping(epoch,track_pair).track1_mean_replay_spikes=mean(remapping(epoch,track_pair).spikes1,2);
            remapping(epoch,track_pair).track2_mean_replay_spikes=mean(remapping(epoch,track_pair).spikes2,2);
           remapping(epoch,track_pair).track1_mean_replay_spikes_nonZero=sum(remapping(epoch,track_pair).spikes1,2)./sum(sign(remapping(epoch,track_pair).spikes1),2); %mean of events with 1 or more spikes
         remapping(epoch,track_pair).track2_mean_replay_spikes_nonZero=sum(remapping(epoch,track_pair).spikes2,2)./sum(sign(remapping(epoch,track_pair).spikes2),2); 
            remapping(epoch,track_pair).track1_mean_replay_rate=mean(remapping(epoch,track_pair).rate1,2);
            remapping(epoch,track_pair).track2_mean_replay_rate=mean(remapping(epoch,track_pair).rate2,2);
            
            index=find(remapping(epoch,track_pair).track1_mean_replay_spikes~=0 & remapping(epoch,track_pair).track2_mean_replay_spikes~=0);  %only analyze neurons that are active during replay.
            
            remapping(epoch,track_pair).fraction_of_cells_not_silent_during_replay=length(index)/length(remapping(epoch,track_pair).track1_mean_replay_spikes);
            remapping(epoch,track_pair).place_field_centre_diff=remapping(epoch,track_pair).place_field_centre_diff(index);
            remapping(epoch,track_pair).track1_mean_replay_rate=remapping(epoch,track_pair).track1_mean_replay_rate(index);
            remapping(epoch,track_pair).track2_mean_replay_rate=remapping(epoch,track_pair).track2_mean_replay_rate(index);
            remapping(epoch,track_pair).track1_mean_replay_spikes=remapping(epoch,track_pair).track1_mean_replay_spikes(index);
            remapping(epoch,track_pair).track2_mean_replay_spikes=remapping(epoch,track_pair).track2_mean_replay_spikes(index);
            remapping(epoch,track_pair).place_fields_BAYESIAN1=remapping(epoch,track_pair).place_fields_BAYESIAN1(index);
            remapping(epoch,track_pair).place_fields_BAYESIAN2=remapping(epoch,track_pair).place_fields_BAYESIAN2(index);
             remapping(epoch,track_pair).track1_mean_replay_spikes_nonZero=remapping(epoch,track_pair).track1_mean_replay_spikes_nonZero(index);
            remapping(epoch,track_pair).track2_mean_replay_spikes_nonZero=remapping(epoch,track_pair).track2_mean_replay_spikes_nonZero(index);
            
            remapping(epoch,track_pair).place_field_diff=(remapping(epoch,track_pair).place_fields_BAYESIAN1-remapping(epoch,track_pair).place_fields_BAYESIAN2)';
            remapping(epoch,track_pair).replay_spike_diff=(remapping(epoch,track_pair).track1_mean_replay_spikes-remapping(epoch,track_pair).track2_mean_replay_spikes);
             remapping(epoch,track_pair).replay_spike_diff_nonZero=remapping(epoch,track_pair).track1_mean_replay_spikes_nonZero-remapping(epoch,track_pair).track2_mean_replay_spikes_nonZero;
           
            remapping(epoch,track_pair).replay_rate_diff=(remapping(epoch,track_pair).track1_mean_replay_rate-remapping(epoch,track_pair).track2_mean_replay_rate);
            
            remapping(epoch,track_pair).tracks_compared=tracks_to_compare;
            
        end
    end
end
end