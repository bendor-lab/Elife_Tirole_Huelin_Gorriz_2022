function  chimeric_track_analysis
load extracted_position.mat
load extracted_replay_events.mat
load extracted_sleep_state.mat
load extracted_clusters.mat

if exist('chimeric')~=7
    mkdir chimeric
end
create_chimeric_track('BAYESIAN');
create_chimeric_track('NORMAL');

%copy relevant data from main folder to chimeric folder
cd chimeric
save extracted_position position
save extracted_replay_events replay reactivations
save extracted_sleep_state sleep_state
save extracted_clusters clusters

batch_analysis('REPLAY');  %run batch_analysis for decoding replay and running shuffles
batch_analysis('ANALYSIS');

cd ..
end


