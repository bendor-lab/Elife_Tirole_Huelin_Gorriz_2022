function  global_remapped_track_analysis(RExp)


load extracted_position.mat
load extracted_replay_events.mat
load extracted_sleep_state.mat
load extracted_clusters.mat
load extracted_place_fields_BAYESIAN.mat
load extracted_place_fields.mat

if exist('global_remapped')~=7
    mkdir global_remapped
end

%copy relevant data from main folder to chimeric folder
cd global_remapped
save extracted_position position
save extracted_replay_events replay reactivations
save extracted_sleep_state sleep_state
save extracted_clusters clusters
save extracted_place_fields_BAYESIAN.mat
save extracted_place_fields.mat

create_global_remapped_track('BAYESIAN');
create_global_remapped_track('NORMAL');

remapping_pipeline('REPLAY',RExp);  %run batch_analysis for decoding replay and running shuffles
remapping_pipeline('ANALYSIS',RExp);
cd ..
end


