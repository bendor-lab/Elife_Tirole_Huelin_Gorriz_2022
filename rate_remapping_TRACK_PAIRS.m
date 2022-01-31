function [remapping, remapping_raw] = rate_remapping_TRACK_PAIRS(folders,varargin)
% Input:
    % folders is [] for individual session, or cell array of folders
    % varargin{1} is for methods, and can be 'spearman' or 'wcorr'
    % varargin{2} as save toggle - can be 1 (default) or 0
% Output: 
    % remapping: difference between variables
    % remapping_raw: variables for each track

if isempty(folders)
%     load('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\folders_to_process_remapping.mat')
    a=pwd;
    folders={a};
end
master_folder= pwd;

if ~isempty(varargin)
    method= varargin{1};
    if length(varargin)>1
        save_option = varargin{2};
    else
        save_option= 1;
    end
    disp(['running for... ' varargin{1}])
else
    method= 'wcorr';
    save_option = 1;
end


% Allocate variables in two structures
track_pair = 1; %not removing but would only be used if we ever have more than two tracks
for epoch = 1 : 7 % PRE & POST, sleep and awake   
    remapping(epoch,track_pair).experiment=[];
    remapping(epoch,track_pair).experiment_all= [];
    remapping(epoch,track_pair).folder= [];
    remapping(epoch,track_pair).epoch= [];
    remapping(epoch,track_pair).common_good_cells=[];
    remapping(epoch,track_pair).ID_active_cells_during_replay=[];
    remapping(epoch,track_pair).new_ID =[];
    remapping(epoch,track_pair).PRE_to_POST_active_cells=[];
    remapping(epoch,track_pair).fraction_of_active_cells_during_replay=[];
    remapping(epoch,track_pair).place_field_diff=[];
    remapping(epoch,track_pair).log_place_field_diff= [];
    remapping(epoch,track_pair).place_field_centre_diff=[];
    remapping(epoch,track_pair).replay_spike_diff=[];
    remapping(epoch,track_pair).replay_spike_diff_nonZero=[];
    remapping(epoch,track_pair).replay_spike_median_diff_nonZero=[];
    remapping(epoch,track_pair).replay_rate_diff=[];
    remapping(epoch,track_pair).log_replay_rate_diff=[];
    remapping(epoch,track_pair).median_replay_rate_diff=[];
    remapping(epoch,track_pair).track_proportion_events_diff=[];
    remapping(epoch,track_pair).track_mean_rate_diff= [];
    remapping(epoch,track_pair).mean_max_FR_replay_diff= [];
    
    remapping_raw(epoch,track_pair).experiment=[];
    remapping_raw(epoch,track_pair).folder= [];
    remapping_raw(epoch,track_pair).epoch= [];
    remapping_raw(epoch,track_pair).common_good_cells=[];
    remapping_raw(epoch,track_pair).ID_active_cells_during_replay=[];
    remapping_raw(epoch,track_pair).new_ID =[];
    remapping_raw(epoch,track_pair).PRE_to_POST_active_cells=[];
    remapping_raw(epoch,track_pair).raw_peak_BAYESIAN_plfield_1=[];
    remapping_raw(epoch,track_pair).raw_peak_BAYESIAN_plfield_2=[];
    remapping_raw(epoch,track_pair).track1_mean_replay_spikes=[];
    remapping_raw(epoch,track_pair).track2_mean_replay_spikes=[];
    remapping_raw(epoch,track_pair).track1_mean_replay_spikes_nonZero=[];
    remapping_raw(epoch,track_pair).track2_mean_replay_spikes_nonZero=[];
    remapping_raw(epoch,track_pair).track1_median_replay_spikes_nonZero=[];
    remapping_raw(epoch,track_pair).track2_median_replay_spikes_nonZero=[];
    remapping_raw(epoch,track_pair).track1_median_replay_rate=[];
    remapping_raw(epoch,track_pair).track2_median_replay_rate=[];
    remapping_raw(epoch,track_pair).replay_spikes1=[];
    remapping_raw(epoch,track_pair).replay_spikes2=[];
    remapping_raw(epoch,track_pair).replay_rate1=[];
    remapping_raw(epoch,track_pair).replay_rate2=[];
    remapping_raw(epoch,track_pair).proportion_events_cell_active_track1= [];
    remapping_raw(epoch,track_pair).proportion_events_cell_active_track2= [];
    remapping_raw(epoch,track_pair).mean_rate_T1= [];
    remapping_raw(epoch,track_pair).mean_rate_T2= [];
    remapping_raw(epoch,track_pair).proportion_events_T1= [];
    remapping_raw(epoch,track_pair).proportion_events_T2= [];
    remapping_raw(epoch,track_pair).max_instantaneous_FR_1= [];
    remapping_raw(epoch,track_pair).max_instantaneous_FR_2= [];
    remapping_raw(epoch,track_pair).track1_mean_replay_inst_FR_nonZero = [];
    remapping_raw(epoch,track_pair).track2_mean_replay_inst_FR_nonZero = [];
    
%       remapping(epoch,track_pair).peak_FR_T1=[];
%       remapping(epoch,track_pair).peak_FR_T2=[];    
end
cd(master_folder);

for i = 1 : length(folders)
    cd(folders{i});
    data = compare_replay_across_tracks(method);
    for this_epoch=1:7
        % Get new cell IDs
        if ~isempty(data(this_epoch).ID_active_cells_during_replay) %if there's PRE session
            data(this_epoch).new_ID = create_new_remapping_IDs(method,folders{i},data(this_epoch).ID_active_cells_during_replay);
        end
    end
%     data(2).new_ID = create_new_remapping_IDs(method,folders{i},data(2).ID_active_cells_during_replay);
    % Find common cells that are active from PRE to POST and save
    [data(1).PRE_to_POST_active_cells, data(2).PRE_to_POST_active_cells ] = deal(intersect(data(2).new_ID,data(1).new_ID));
    [data(3).PRE_to_POST_active_cells, data(4).PRE_to_POST_active_cells ] = deal(intersect(data(3).new_ID,data(4).new_ID));
    [data(5).PRE_to_POST_active_cells, data(6).PRE_to_POST_active_cells ] = deal(intersect(data(5).new_ID,data(6).new_ID));

    for epoch = 1 : size(data,1)
        for track_pair = 1 : size(data,2)
            remapping(epoch,track_pair).folder = [remapping(epoch,track_pair).folder folders(i)];
            remapping_raw(epoch,track_pair).folder = [remapping_raw(epoch,track_pair).folder folders(i)];
            if epoch==1
                remapping(epoch,track_pair).epoch= {'sleepPRE'};
                remapping_raw(epoch,track_pair).epoch= {'sleepPRE'};                
            elseif epoch ==2
                remapping(epoch,track_pair).epoch= {'sleepPOST'};
                remapping_raw(epoch,track_pair).epoch= {'sleepPOST'};
            elseif epoch ==3
                remapping(epoch,track_pair).epoch= {'awakePRE'};
                remapping_raw(epoch,track_pair).epoch= {'awakePRE'};
            elseif epoch ==4
                remapping(epoch,track_pair).epoch= {'awakePOST'};
                remapping_raw(epoch,track_pair).epoch= {'awakePOST'};
            elseif epoch==5
                remapping(epoch,track_pair).epoch= {'PRE'};
                remapping_raw(epoch,track_pair).epoch= {'PRE'};
            elseif epoch==6
                remapping(epoch,track_pair).epoch= {'POST'};
                remapping_raw(epoch,track_pair).epoch= {'POST'};
            elseif epoch==7
                remapping(epoch,track_pair).epoch= {'RUN'};
                remapping_raw(epoch,track_pair).epoch= {'RUN'};
            end
            remapping(epoch,track_pair).experiment = [remapping(epoch,track_pair).experiment; (i*ones(size(data(epoch,track_pair).place_field_diff)))];
            remapping_raw(epoch,track_pair).experiment = [remapping_raw(epoch,track_pair).experiment; (i*ones(size(data(epoch,track_pair).place_field_diff)))];
            remapping(epoch,track_pair).experiment_all = [remapping(epoch,track_pair).experiment_all; (i*ones(size(data(epoch,track_pair).good_cells')))];
            remapping(epoch,track_pair).common_good_cells = [remapping(epoch,track_pair).common_good_cells; data(epoch,track_pair).good_cells'];
            remapping_raw(epoch,track_pair).common_good_cells = [remapping_raw(epoch,track_pair).common_good_cells; data(epoch,track_pair).good_cells'];
            remapping(epoch,track_pair).ID_active_cells_during_replay = [remapping(epoch,track_pair).ID_active_cells_during_replay; data(epoch,track_pair).ID_active_cells_during_replay'];
            remapping_raw(epoch,track_pair).ID_active_cells_during_replay = [remapping_raw(epoch,track_pair).ID_active_cells_during_replay; data(epoch,track_pair).ID_active_cells_during_replay'];
            remapping(epoch,track_pair).new_ID = [remapping(epoch,track_pair).new_ID data(epoch,track_pair).new_ID];
            remapping_raw(epoch,track_pair).new_ID = [remapping_raw(epoch,track_pair).new_ID data(epoch,track_pair).new_ID];
            remapping(epoch,track_pair).fraction_of_active_cells_during_replay = [remapping(epoch,track_pair).fraction_of_active_cells_during_replay; data(epoch,track_pair).fraction_of_active_cells_during_replay];
            remapping(epoch,track_pair).PRE_to_POST_active_cells = [remapping(epoch,track_pair).PRE_to_POST_active_cells; data(epoch,track_pair).PRE_to_POST_active_cells'];
            remapping_raw(epoch,track_pair).PRE_to_POST_active_cells = [remapping_raw(epoch,track_pair).PRE_to_POST_active_cells; data(epoch,track_pair).PRE_to_POST_active_cells'];

            % Remapping- diff tracks structure
            remapping(epoch,track_pair).place_field_diff =[remapping(epoch,track_pair).place_field_diff; data(epoch,track_pair).place_field_diff];
            remapping(epoch,track_pair).log_place_field_diff = [remapping(epoch,track_pair).log_place_field_diff; data(epoch,track_pair).place_field_diff];
            remapping(epoch,track_pair).place_field_centre_diff = [remapping(epoch,track_pair).place_field_centre_diff; data(epoch,track_pair).place_field_centre_diff];
            remapping(epoch,track_pair).replay_spike_diff = [remapping(epoch,track_pair).replay_spike_diff; data(epoch,track_pair).replay_spike_diff];
            remapping(epoch,track_pair).replay_spike_diff_nonZero = [remapping(epoch,track_pair).replay_spike_diff_nonZero; data(epoch,track_pair).replay_spike_diff_nonZero];
            remapping(epoch,track_pair).replay_spike_median_diff_nonZero = [remapping(epoch,track_pair).replay_spike_median_diff_nonZero; data(epoch,track_pair).replay_spike_median_diff_nonZero];
            remapping(epoch,track_pair).median_replay_rate_diff = [remapping(epoch,track_pair).median_replay_rate_diff; data(epoch,track_pair).median_replay_rate_diff];
            remapping(epoch,track_pair).replay_rate_diff = [remapping(epoch,track_pair).replay_rate_diff; data(epoch,track_pair).replay_rate_diff];
            remapping(epoch,track_pair).log_replay_rate_diff = [remapping(epoch,track_pair).log_replay_rate_diff; data(epoch,track_pair).log_replay_rate_diff];
            remapping(epoch,track_pair).track_proportion_events_diff = [remapping(epoch,track_pair).track_proportion_events_diff data(epoch,track_pair).proportion_events_diff];
            remapping(epoch,track_pair).track_mean_rate_diff = [remapping(epoch,track_pair).track_mean_rate_diff; data(epoch,track_pair).mean_rate_diff'];
            remapping(epoch,track_pair).mean_max_FR_replay_diff= [remapping(epoch,track_pair).mean_max_FR_replay_diff;  data(epoch,track_pair).mean_max_FR_replay_diff];
            
            % Remapping - track data structure
            remapping_raw(epoch,track_pair).raw_peak_BAYESIAN_plfield_1 = [remapping_raw(epoch,track_pair).raw_peak_BAYESIAN_plfield_1 data(epoch,track_pair).raw_peak_BAYESIAN_plfield_1];
            remapping_raw(epoch,track_pair).raw_peak_BAYESIAN_plfield_2 = [remapping_raw(epoch,track_pair).raw_peak_BAYESIAN_plfield_2 data(epoch,track_pair).raw_peak_BAYESIAN_plfield_2];
            remapping_raw(epoch,track_pair).track1_mean_replay_spikes = [remapping_raw(epoch,track_pair).track1_mean_replay_spikes; data(epoch,track_pair).track1_mean_replay_spikes];
            remapping_raw(epoch,track_pair).track2_mean_replay_spikes = [remapping_raw(epoch,track_pair).track2_mean_replay_spikes; data(epoch,track_pair).track2_mean_replay_spikes];
            remapping_raw(epoch,track_pair).track1_median_replay_rate = [remapping_raw(epoch,track_pair).track1_median_replay_rate; data(epoch,track_pair).track1_median_replay_rate];
            remapping_raw(epoch,track_pair).track2_median_replay_rate=  [remapping_raw(epoch,track_pair).track2_median_replay_rate; data(epoch,track_pair).track2_median_replay_rate];
            remapping_raw(epoch,track_pair).track1_mean_replay_spikes_nonZero = [remapping_raw(epoch,track_pair).track1_mean_replay_spikes_nonZero; data(epoch,track_pair).track1_mean_replay_spikes_nonZero];
            remapping_raw(epoch,track_pair).track2_mean_replay_spikes_nonZero = [remapping_raw(epoch,track_pair).track2_mean_replay_spikes_nonZero; data(epoch,track_pair).track2_mean_replay_spikes_nonZero];
            remapping_raw(epoch,track_pair).track1_median_replay_spikes_nonZero = [remapping_raw(epoch,track_pair).track1_median_replay_spikes_nonZero; data(epoch,track_pair).track1_median_replay_spikes_nonZero];
            remapping_raw(epoch,track_pair).track2_median_replay_spikes_nonZero = [remapping_raw(epoch,track_pair).track2_median_replay_spikes_nonZero; data(epoch,track_pair).track2_median_replay_spikes_nonZero];
            remapping_raw(epoch,track_pair).proportion_events_cell_active_track1 = [remapping_raw(epoch,track_pair).proportion_events_cell_active_track1; data(epoch,track_pair).proportion_events_cell_active_track1];
            remapping_raw(epoch,track_pair).proportion_events_cell_active_track2 = [remapping_raw(epoch,track_pair).proportion_events_cell_active_track2; data(epoch,track_pair).proportion_events_cell_active_track2];
            remapping_raw(epoch,track_pair).mean_rate_T1 = [remapping_raw(epoch,track_pair).mean_rate_T1; data(epoch,track_pair).mean_rate_T1'];
            remapping_raw(epoch,track_pair).mean_rate_T2 = [remapping_raw(epoch,track_pair).mean_rate_T2; data(epoch,track_pair).mean_rate_T2'];
            remapping_raw(epoch,track_pair).proportion_events_T1 = [remapping_raw(epoch,track_pair).proportion_events_T1 data(epoch,track_pair).proportion_events_T1];
            remapping_raw(epoch,track_pair).proportion_events_T2 = [remapping_raw(epoch,track_pair).proportion_events_T2 data(epoch,track_pair).proportion_events_T2];
            remapping_raw(epoch,track_pair).replay_spikes1 = [remapping_raw(epoch,track_pair).replay_spikes1 {data(epoch,track_pair).spikes1}];
            remapping_raw(epoch,track_pair).replay_spikes2 = [remapping_raw(epoch,track_pair).replay_spikes2 {data(epoch,track_pair).spikes2}];
            remapping_raw(epoch,track_pair).replay_rate1 = [remapping_raw(epoch,track_pair).replay_rate1 {data(epoch,track_pair).rate1}];
            remapping_raw(epoch,track_pair).replay_rate2 = [remapping_raw(epoch,track_pair).replay_rate2 {data(epoch,track_pair).rate2}];
            remapping_raw(epoch,track_pair).max_instantaneous_FR_1 = [remapping_raw(epoch,track_pair).max_instantaneous_FR_1 {data(epoch,track_pair).max_instantaneous_FR_1}];
            remapping_raw(epoch,track_pair).max_instantaneous_FR_2 = [remapping_raw(epoch,track_pair).max_instantaneous_FR_2 {data(epoch,track_pair).max_instantaneous_FR_2}];
            remapping_raw(epoch,track_pair).track1_mean_replay_inst_FR_nonZero = [remapping_raw(epoch,track_pair).track1_mean_replay_inst_FR_nonZero; data(epoch,track_pair).track1_mean_replay_inst_FR_nonZero];
            remapping_raw(epoch,track_pair).track2_mean_replay_inst_FR_nonZero = [remapping_raw(epoch,track_pair).track2_mean_replay_inst_FR_nonZero; data(epoch,track_pair).track2_mean_replay_inst_FR_nonZero];

%           remapping(epoch,track_pair).peak_FR_T1=  [remapping(epoch,track_pair).peak_FR_T1; data(epoch,track_pair).peak_FR_T1'];
%           remapping(epoch,track_pair).peak_FR_T2=  [remapping(epoch,track_pair).peak_FR_T2; data(epoch,track_pair).peak_FR_T2'];
             
        end
    end
   cd(master_folder);
end

if ~isempty(folders) && save_option == 1
    switch method
        case 'wcorr'
            save rate_remapping_analysis_TRACK_PAIRS_wcorr remapping remapping_raw
        case 'spearman'
            save rate_remapping_analysis_TRACK_PAIRS_spearman remapping remapping_raw
        case 'control_fixed_spike'
            save rate_remapping_analysis_TRACK_PAIRS_wcorr_FIXED remapping remapping_raw
        case 'rate_detection_control'
            save rate_remapping_analysis_TRACK_PAIRS_wcorr_RATE remapping remapping_raw
        case 'rate_intrinsic_bias_control'
            save rate_remapping_analysis_TRACK_PAIRS_wcorr_INTRINSIC_RATE remapping remapping_raw
        case 'replay_rate_shuffle_control'
            save rate_remapping_analysis_TRACK_PAIRS_wcorr_REPLAY_RATE remapping remapping_raw
        case 'replay_rate_shuffle_detection'
            save rate_remapping_analysis_TRACK_PAIRS_wcorr_REPLAY_RATE_DETECTION remapping remapping_raw
        otherwise
            save rate_remapping_analysis_TRACK_PAIRS remapping remapping_raw
    end
end

end


function remapping= compare_replay_across_tracks(method)

switch method
    case 'wcorr'
        if exist('significant_replay_events_wcorr.mat')==2
            load('significant_replay_events_wcorr.mat');
        elseif exist('significant_replay_events_wcorr_individual_exposures.mat')==2
            load('significant_replay_events_wcorr_individual_exposures');
        end
        load('sorted_replay_wcorr');
    case 'spearman'
        if exist('significant_replay_events_spearman.mat')==2
            load('significant_replay_events_spearman.mat');
        elseif exist('significant_replay_events_spearman_individual_exposures.mat')==2
            load('significant_replay_events_spearman_individual_exposures.mat');
        end
        load('sorted_replay_spearman');
    case 'control_fixed_spike'
        if exist('significant_replay_events_wcorr_FIXED.mat')==2
            load('significant_replay_events_wcorr_FIXED.mat');
        elseif exist('significant_replay_events_wcorr_individual_exposures_FIXED.mat')==2
            load('significant_replay_events_wcorr_individual_exposures_FIXED');
        end
        load('sorted_replay_wcorr_FIXED');
     case 'rate_detection_control'
        if exist('significant_replay_events_wcorr_RATE.mat')==2
            load('significant_replay_events_wcorr_RATE.mat');
        elseif exist('significant_replay_events_wcorr_individual_exposures_RATE.mat')==2
            load('significant_replay_events_wcorr_individual_exposures_RATE');
        end
        load('sorted_replay_wcorr_RATE');
        data_folder = strfind(pwd,'\');
        fname = pwd; fname = fname(data_folder(end)+1:end);
        load(['..\..\..\' fname '\extracted_place_fields_BAYESIAN.mat']); % load real place fields       
    case 'rate_intrinsic_bias_control'
        data_folder = strfind(pwd,'\'); %load real events
        fname = pwd; fname = fname(data_folder(end)+1:end);
        if contains(fname,'TEMP') % when running shuffle distribution
           fname = fname(6:end); 
        end
        if exist(['..\..\..\' fname '\significant_replay_events_wcorr.mat'])==2
            load(['..\..\..\' fname '\significant_replay_events_wcorr.mat']);
        elseif exist(['..\..\..\' fname '\significant_replay_events_wcorr_individual_exposures.mat'])==2
            load(['..\..\..\' fname '\significant_replay_events_wcorr_individual_exposures']);
        end
        load(['..\..\..\' fname '\sorted_replay_wcorr.mat']);   
    case 'replay_rate_shuffle_control'
        data_folder = strfind(pwd,'\');
        fname = pwd; fname = fname(data_folder(end)+1:end);
        if contains(fname,'TEMP_') % when running shuffle distribution
           fname = fname(6:end); 
        elseif contains(fname,'TEMP2_')
            fname = fname(7:end); 
        end
        if exist(['..\..\..\' fname '\significant_replay_events_wcorr.mat'])==2 % load original events
            load(['..\..\..\' fname '\significant_replay_events_wcorr.mat']);
        elseif exist(['..\..\..\' fname '\significant_replay_events_wcorr_individual_exposures.mat'])==2
            load(['..\..\..\' fname '\significant_replay_events_wcorr_individual_exposures']);
        end
        load(['..\..\..\' fname '\sorted_replay_wcorr.mat']);
        load('replayEvents_bayesian_spike_count.mat'); % to get scaling factors for spikes
   case 'replay_rate_shuffle_detection'
        data_folder = strfind(pwd,'\');
        fname = pwd; fname = fname(data_folder(end)+1:end);
        if exist(['..\..\..\' fname '\significant_replay_events_wcorr.mat'])==2 % load original spikes
            load(['..\..\..\' fname '\significant_replay_events_wcorr.mat']);
        elseif exist(['..\..\..\' fname '\significant_replay_events_wcorr_individual_exposures.mat'])==2
            load(['..\..\..\' fname '\significant_replay_events_wcorr_individual_exposures']);
        end
        load('sorted_replay_wcorr'); % but load new events
     otherwise
        load('significant_replay_events');
        load('sorted_replay');
end

if ~exist('place_fields_BAYESIAN','var')
    load('extracted_place_fields_BAYESIAN');
end 

number_of_tracks=length(sorted_replay);
track_combinations=[1 2];
% if number_of_tracks==2
%     track_combinations=[1 2];
% elseif number_of_tracks==3
%     track_combinations=[1 2; 2 3; 3 1];
% elseif number_of_tracks==4
%     track_combinations=[1 2; 2 3; 3 4; 4 1; 1 3; 2 4];
% end

% create gaussian kernel
n_sigma= 0.1; % 100ms
gauss_ker= gausswin(n_sigma/(1/1000));
% gauss_ker= gauss_ker./sum(gauss_ker); % this would ideally be removed, bad normalisation

for track_pair = 1:size(track_combinations,1)
    tracks_to_compare=track_combinations(track_pair,:);
    for epoch = 1:7  %sleep PRE &  POST, awake PRE & POST, PRE & POST, RUN
        if epoch==1
            index1  = sorted_replay(tracks_to_compare(1)).index.sleepPRE; %indices replay events PRE for T1
            events1 = significant_replay_events.track(tracks_to_compare(1)); %info from events in T1
            index2  = sorted_replay(tracks_to_compare(2)).index.sleepPRE; %indices replay events PRE for T2
            events2 = significant_replay_events.track(tracks_to_compare(2)); %info from events in T2
            ref_index1= sorted_replay(tracks_to_compare(1)).ref_index.sleepPRE;
            ref_index2= sorted_replay(tracks_to_compare(2)).ref_index.sleepPRE;
        elseif epoch==2
            index1  = sorted_replay(tracks_to_compare(1)).index.sleepPOST; %indices replay events POST for T1
            events1 = significant_replay_events.track(tracks_to_compare(1));
            index2  = sorted_replay(tracks_to_compare(2)).index.sleepPOST; %indices replay events POST for T2
            events2 = significant_replay_events.track(tracks_to_compare(2));
            ref_index1= sorted_replay(tracks_to_compare(1)).ref_index.sleepPOST;
            ref_index2= sorted_replay(tracks_to_compare(2)).ref_index.sleepPOST;
        elseif epoch==3
            index1  = sorted_replay(tracks_to_compare(1)).index.awakePRE; 
            events1 = significant_replay_events.track(tracks_to_compare(1)); 
            index2  = sorted_replay(tracks_to_compare(2)).index.awakePRE; 
            events2 = significant_replay_events.track(tracks_to_compare(2)); 
            ref_index1= sorted_replay(tracks_to_compare(1)).ref_index.awakePRE;
            ref_index2= sorted_replay(tracks_to_compare(2)).ref_index.awakePRE;
        elseif epoch==4
            index1  = sorted_replay(tracks_to_compare(1)).index.awakePOST; 
            events1 = significant_replay_events.track(tracks_to_compare(1)); 
            index2  = sorted_replay(tracks_to_compare(2)).index.awakePOST; 
            events2 = significant_replay_events.track(tracks_to_compare(2)); 
            ref_index1= sorted_replay(tracks_to_compare(1)).ref_index.awakePOST;
            ref_index2= sorted_replay(tracks_to_compare(2)).ref_index.awakePOST;
        elseif epoch==5
            index1  = sorted_replay(tracks_to_compare(1)).index.PRE; 
            events1 = significant_replay_events.track(tracks_to_compare(1)); 
            index2  = sorted_replay(tracks_to_compare(2)).index.PRE; 
            events2 = significant_replay_events.track(tracks_to_compare(2)); 
            ref_index1= sorted_replay(tracks_to_compare(1)).ref_index.PRE;
            ref_index2= sorted_replay(tracks_to_compare(2)).ref_index.PRE;
        elseif epoch==6
            index1  = sorted_replay(tracks_to_compare(1)).index.POST; 
            events1 = significant_replay_events.track(tracks_to_compare(1)); 
            index2  = sorted_replay(tracks_to_compare(2)).index.POST; 
            events2 = significant_replay_events.track(tracks_to_compare(2)); 
            ref_index1= sorted_replay(tracks_to_compare(1)).ref_index.POST;
            ref_index2= sorted_replay(tracks_to_compare(2)).ref_index.POST;
        elseif epoch == 7
            index1  = sorted_replay(tracks_to_compare(1)).index.track(tracks_to_compare(1)).behaviour; 
            events1 = significant_replay_events.track(tracks_to_compare(1)); 
            index2  = sorted_replay(tracks_to_compare(2)).index.track(tracks_to_compare(2)).behaviour; 
            events2 = significant_replay_events.track(tracks_to_compare(2)); 
            ref_index1= sorted_replay(tracks_to_compare(1)).ref_index.track(tracks_to_compare(1)).behaviour;
            ref_index2= sorted_replay(tracks_to_compare(2)).ref_index.track(tracks_to_compare(2)).behaviour;
        end
        
        remapping(epoch,track_pair).fraction_of_active_cells_during_replay =[];
        remapping(epoch,track_pair).place_field_diff=[];
        remapping(epoch,track_pair).replay_spike_diff=[];
        remapping(epoch,track_pair).replay_rate_diff=[];
        remapping(epoch,track_pair).place_field_centre_diff=[];
        remapping(epoch,track_pair).tracks_compared=[];
        remapping(epoch,track_pair).track1_mean_replay_spikes=[];
        remapping(epoch,track_pair).track2_mean_replay_spikes=[];
        remapping(epoch,track_pair).track1_median_replay_rate=[];
        remapping(epoch,track_pair).track2_median_replay_rate=[];
        remapping(epoch,track_pair).median_replay_rate_diff=[];
        remapping(epoch,track_pair).raw_peak_BAYESIAN_plfield_1=[];
        remapping(epoch,track_pair).raw_peak_BAYESIAN_plfield_2=[]; 
        
        % Find BAYESIAN good cells in common for both tracks
        remapping(epoch,track_pair).good_cells= intersect(place_fields_BAYESIAN.track(tracks_to_compare(1)).good_cells, place_fields_BAYESIAN.track(tracks_to_compare(2)).good_cells);
        
         % Calculate number of spikes in replay and divide by duration (based on 20ms binning, not exact time due to binning of spikes)
        if length(index1)>0 & length(index2)>0 %if there's events for both tracks 
            for j=1:length(index1) %for each event in T1
                for i=1:length(remapping(epoch,track_pair).good_cells) %for each common cell
                    % Find # spikes from each common cell in this event 
                    remapping(epoch,track_pair).spikes1(i,j) = length(find(events1.spikes{index1(j)}(:,1)==remapping(epoch,track_pair).good_cells(i)));
                    % Find FR from each common cell in this event
                    remapping(epoch,track_pair).rate1(i,j) = length(find(events1.spikes{index1(j)}(:,1)==remapping(epoch,track_pair).good_cells(i)))/events1.event_duration(index1(j));
                    % Find max instantaneous firing rate during replay
                    % create edges with 1ms bins, add some length of filter before and after for filtering ease
                    replay_time= events1.event_times(index1(j));
                    replay_dur= events1.event_duration(index1(j));
                    replay_spikes= events1.spikes{index1(j)};
                    ts_event_edges= min(replay_time-(replay_dur/2),min(replay_spikes(:,2)))-n_sigma :1/1000: max(replay_time+(replay_dur/2),max(replay_spikes(:,2)))+n_sigma;
                    binned_spike_train= histcounts(events1.spikes{index1(j)}(events1.spikes{index1(j)}(:,1)==remapping(epoch,track_pair).good_cells(i),2),ts_event_edges);
                    if strcmp(method,'replay_rate_shuffle_control')
                        cell_idx= place_fields_BAYESIAN.good_place_cells==remapping(epoch,track_pair).good_cells(i);
                        event_ref_idx= ref_index1(j);
%                         rate_event= length(replay_spikes(events1.spikes{index1(j)}(:,1)==remapping(epoch,track_pair).good_cells(i)))./replay_dur;
                        binned_spike_train= binned_spike_train.*(replayEvents_bayesian_spike_count.scaling_factors(cell_idx,event_ref_idx));
                    end
                    instantaneous_FR = filter(gauss_ker,1,binned_spike_train);
                    remapping(epoch,track_pair).max_instantaneous_FR_1(i,j) = max(instantaneous_FR); % we only care about amplitude, not timing
                end
            end
            for j=1:length(index2) %for each event in T2
                for i=1:length(remapping(epoch,track_pair).good_cells)
                    remapping(epoch,track_pair).spikes2(i,j)=length(find(events2.spikes{index2(j)}(:,1)==remapping(epoch,track_pair).good_cells(i)));
                    remapping(epoch,track_pair).rate2(i,j)=length(find(events2.spikes{index2(j)}(:,1)==remapping(epoch,track_pair).good_cells(i)))/events2.event_duration(index2(j));
                    replay_time= events2.event_times(index2(j));
                    replay_dur= events2.event_duration(index2(j));
                    replay_spikes= events2.spikes{index2(j)};
                    ts_event_edges= min(replay_time-(replay_dur/2),min(replay_spikes(:,2)))-n_sigma :1/1000: max(replay_time+(replay_dur/2),max(replay_spikes(:,2)))+n_sigma;
                    binned_spike_train= histcounts(events2.spikes{index2(j)}(events2.spikes{index2(j)}(:,1)==remapping(epoch,track_pair).good_cells(i),2),ts_event_edges);
                    if strcmp(method,'replay_rate_shuffle_control')
                        cell_idx= place_fields_BAYESIAN.good_place_cells==remapping(epoch,track_pair).good_cells(i);
                        event_ref_idx= ref_index2(j);
                        binned_spike_train= binned_spike_train.*(replayEvents_bayesian_spike_count.scaling_factors(cell_idx,event_ref_idx));
                    end
                    instantaneous_FR = filter(gauss_ker,1,binned_spike_train);
                    remapping(epoch,track_pair).max_instantaneous_FR_2(i,j) = max(instantaneous_FR); 
                end
            end
            
            %%%% Track 1 - track 2 comparison of only good place cells (on both tracks)
            
            % Get raw peak FR
            remapping(epoch,track_pair).raw_peak_BAYESIAN_plfield_1= place_fields_BAYESIAN.track(tracks_to_compare(1)).raw_peak(remapping(epoch,track_pair).good_cells);
            remapping(epoch,track_pair).raw_peak_BAYESIAN_plfield_2= place_fields_BAYESIAN.track(tracks_to_compare(2)).raw_peak(remapping(epoch,track_pair).good_cells);
            % get mean rate on track
            remapping(epoch,track_pair).mean_rate_T1= place_fields_BAYESIAN.track(tracks_to_compare(1)).mean_rate_track(remapping(epoch,track_pair).good_cells);
            remapping(epoch,track_pair).mean_rate_T2= place_fields_BAYESIAN.track(tracks_to_compare(2)).mean_rate_track(remapping(epoch,track_pair).good_cells);
            remapping(epoch,track_pair).mean_rate_diff= remapping(epoch,track_pair).mean_rate_T1- remapping(epoch,track_pair).mean_rate_T2;

            % Calculate centre of mass difference 
            remapping(epoch,track_pair).place_field_centre_diff=(abs(place_fields_BAYESIAN.track(tracks_to_compare(1)).centre_of_mass(remapping(epoch,track_pair).good_cells)-...
                place_fields_BAYESIAN.track(tracks_to_compare(2)).centre_of_mass(remapping(epoch,track_pair).good_cells)))';
            % Mean number of spikes across all events for each cell 
            remapping(epoch,track_pair).track1_mean_replay_spikes= mean(remapping(epoch,track_pair).spikes1,2);
            remapping(epoch,track_pair).track2_mean_replay_spikes= mean(remapping(epoch,track_pair).spikes2,2);
            % proportion of events cell is active in
            remapping(epoch,track_pair).proportion_events_cell_active_track1= sum(sign(remapping(epoch,track_pair).spikes1),2)/size(remapping(epoch,track_pair).spikes1,2);
            remapping(epoch,track_pair).proportion_events_cell_active_track2= sum(sign(remapping(epoch,track_pair).spikes2),2)/size(remapping(epoch,track_pair).spikes2,2);   
            % proportion events track1 vs track2
            remapping(epoch,track_pair).proportion_events_T1= size(remapping(epoch,track_pair).spikes1,2)/(size(remapping(epoch,track_pair).spikes1,2) + size(remapping(epoch,track_pair).spikes2,2));
            remapping(epoch,track_pair).proportion_events_T2= size(remapping(epoch,track_pair).spikes2,2)/(size(remapping(epoch,track_pair).spikes1,2) + size(remapping(epoch,track_pair).spikes2,2));
            % Sum of spikes per cell across all events divided by number of events where cell is active
            remapping(epoch,track_pair).track1_mean_replay_spikes_nonZero= sum(remapping(epoch,track_pair).spikes1,2)./sum(sign(remapping(epoch,track_pair).spikes1),2); %mean of events with 1 or more spikes
            remapping(epoch,track_pair).track2_mean_replay_spikes_nonZero= sum(remapping(epoch,track_pair).spikes2,2)./sum(sign(remapping(epoch,track_pair).spikes2),2); 
            % median number of spikes when cell is active
            track1_nonzero_spikes= remapping(epoch,track_pair).spikes1; track1_nonzero_spikes(track1_nonzero_spikes==0)= NaN;
            track2_nonzero_spikes= remapping(epoch,track_pair).spikes2; track2_nonzero_spikes(track2_nonzero_spikes==0)= NaN;
            remapping(epoch,track_pair).track1_median_replay_spikes_nonZero= nanmedian(track1_nonzero_spikes,2);
            remapping(epoch,track_pair).track2_median_replay_spikes_nonZero= nanmedian(track2_nonzero_spikes,2);
            % now replace NaNs with 0
            remapping(epoch,track_pair).track1_median_replay_spikes_nonZero(isnan(remapping(epoch,track_pair).track1_median_replay_spikes_nonZero))= 0;
            remapping(epoch,track_pair).track2_median_replay_spikes_nonZero(isnan(remapping(epoch,track_pair).track2_median_replay_spikes_nonZero))= 0; 
            
            % median max instantaneous FR replay when cell is active
            track1_nonzero_inst_FR= remapping(epoch,track_pair).max_instantaneous_FR_1; track1_nonzero_inst_FR(track1_nonzero_inst_FR==0)= NaN; % replace 0 by NaN
            track2_nonzero_inst_FR= remapping(epoch,track_pair).max_instantaneous_FR_2; track2_nonzero_inst_FR(track2_nonzero_inst_FR==0)= NaN;
            remapping(epoch,track_pair).track1_mean_replay_inst_FR_nonZero= nanmean(track1_nonzero_inst_FR,2);
            remapping(epoch,track_pair).track2_mean_replay_inst_FR_nonZero= nanmean(track2_nonzero_inst_FR,2);
            % now replace NaNs with 0
            remapping(epoch,track_pair).track1_mean_replay_inst_FR_nonZero(isnan(remapping(epoch,track_pair).track1_mean_replay_inst_FR_nonZero))= 0;
            remapping(epoch,track_pair).track2_mean_replay_inst_FR_nonZero(isnan(remapping(epoch,track_pair).track2_mean_replay_inst_FR_nonZero))= 0; 
            % get diff
            remapping(epoch,track_pair).mean_max_FR_replay_diff= remapping(epoch,track_pair).track1_mean_replay_inst_FR_nonZero- remapping(epoch,track_pair).track2_mean_replay_inst_FR_nonZero;
      
            % Median FR per cell across events with FR>0
            % make 0 NaNs so get true median
            track1_nonzero_rates= remapping(epoch,track_pair).rate1; track1_nonzero_rates(track1_nonzero_rates==0)= NaN;
            track2_nonzero_rates= remapping(epoch,track_pair).rate2; track2_nonzero_rates(track2_nonzero_rates==0)= NaN;
            remapping(epoch,track_pair).track1_median_replay_rate= nanmedian(track1_nonzero_rates,2);
            remapping(epoch,track_pair).track2_median_replay_rate= nanmedian(track2_nonzero_rates,2);
            % now replace NaNs zith 0
            remapping(epoch,track_pair).track1_median_replay_rate(isnan(remapping(epoch,track_pair).track1_median_replay_rate))= 0;
            remapping(epoch,track_pair).track2_median_replay_rate(isnan(remapping(epoch,track_pair).track2_median_replay_rate))= 0;
            % get difference
            remapping(epoch,track_pair).median_replay_rate_diff= remapping(epoch,track_pair).track1_median_replay_rate - remapping(epoch,track_pair).track2_median_replay_rate;
            
            % Find cells that have a mean FR >0 during replay in both T1 and T2 & save cell IDs
            index = find(remapping(epoch,track_pair).track1_mean_replay_spikes~=0 & remapping(epoch,track_pair).track2_mean_replay_spikes~=0);  %only analyze neurons that are active during replay.
            remapping(epoch,track_pair).ID_active_cells_during_replay = remapping(epoch,track_pair).good_cells(index);
     
            % Fraction of cells from common good cells that are active in replay
            remapping(epoch,track_pair).fraction_of_active_cells_during_replay = length(index)/length(remapping(epoch,track_pair).track1_mean_replay_spikes);
            % Remove silent cells
            remapping(epoch,track_pair).place_field_centre_diff=remapping(epoch,track_pair).place_field_centre_diff(index);
            remapping(epoch,track_pair).track1_median_replay_rate=remapping(epoch,track_pair).track1_median_replay_rate(index);
            remapping(epoch,track_pair).track2_median_replay_rate=remapping(epoch,track_pair).track2_median_replay_rate(index);
            remapping(epoch,track_pair).track1_mean_replay_spikes=remapping(epoch,track_pair).track1_mean_replay_spikes(index);
            remapping(epoch,track_pair).track2_mean_replay_spikes=remapping(epoch,track_pair).track2_mean_replay_spikes(index);
            remapping(epoch,track_pair).raw_peak_BAYESIAN_plfield_1=remapping(epoch,track_pair).raw_peak_BAYESIAN_plfield_1(index);
            remapping(epoch,track_pair).raw_peak_BAYESIAN_plfield_2=remapping(epoch,track_pair).raw_peak_BAYESIAN_plfield_2(index);
            remapping(epoch,track_pair).mean_rate_T1= remapping(epoch,track_pair).mean_rate_T1(index);
            remapping(epoch,track_pair).mean_rate_T2= remapping(epoch,track_pair).mean_rate_T2(index);
            remapping(epoch,track_pair).mean_rate_diff= remapping(epoch,track_pair).mean_rate_diff(index);
            remapping(epoch,track_pair).track1_mean_replay_spikes_nonZero=remapping(epoch,track_pair).track1_mean_replay_spikes_nonZero(index);
            remapping(epoch,track_pair).track2_mean_replay_spikes_nonZero=remapping(epoch,track_pair).track2_mean_replay_spikes_nonZero(index);
            remapping(epoch,track_pair).track1_median_replay_spikes_nonZero=remapping(epoch,track_pair).track1_median_replay_spikes_nonZero(index);
            remapping(epoch,track_pair).track2_median_replay_spikes_nonZero=remapping(epoch,track_pair).track2_median_replay_spikes_nonZero(index);
            remapping(epoch,track_pair).proportion_events_cell_active_track1= remapping(epoch,track_pair).proportion_events_cell_active_track1(index);
            remapping(epoch,track_pair).proportion_events_cell_active_track2= remapping(epoch,track_pair).proportion_events_cell_active_track2(index);
            remapping(epoch,track_pair).max_instantaneous_FR_1= remapping(epoch,track_pair).max_instantaneous_FR_1(index,:);
            remapping(epoch,track_pair).max_instantaneous_FR_2= remapping(epoch,track_pair).max_instantaneous_FR_2(index,:);
            remapping(epoch,track_pair).track1_mean_replay_inst_FR_nonZero = remapping(epoch,track_pair).track1_mean_replay_inst_FR_nonZero(index,:);
            remapping(epoch,track_pair).track2_mean_replay_inst_FR_nonZero = remapping(epoch,track_pair).track2_mean_replay_inst_FR_nonZero(index,:);

            % Difference between T1-T2 raw peak for each cell during RUN
            remapping(epoch,track_pair).place_field_diff=(remapping(epoch,track_pair).raw_peak_BAYESIAN_plfield_1-remapping(epoch,track_pair).raw_peak_BAYESIAN_plfield_2)';
            remapping(epoch,track_pair).log_place_field_diff=(log(remapping(epoch,track_pair).raw_peak_BAYESIAN_plfield_1)-log(remapping(epoch,track_pair).raw_peak_BAYESIAN_plfield_2))';

            % Difference between T1-T2 mean number of spikes during REPLAY
            remapping(epoch,track_pair).replay_spike_diff=(remapping(epoch,track_pair).track1_mean_replay_spikes-remapping(epoch,track_pair).track2_mean_replay_spikes);
            % Difference between T1-T2 mean number of spikes during REPLAY events where cell is active
            remapping(epoch,track_pair).replay_spike_diff_nonZero=remapping(epoch,track_pair).track1_mean_replay_spikes_nonZero-remapping(epoch,track_pair).track2_mean_replay_spikes_nonZero;    
            % Difference between T1-T2 median number of spikes during REPLAY events where cell is active
            remapping(epoch,track_pair).replay_spike_median_diff_nonZero=remapping(epoch,track_pair).track1_median_replay_spikes_nonZero-remapping(epoch,track_pair).track2_median_replay_spikes_nonZero;    
            % Difference median replay rate T1-T2 during REPLAY events where cell is active
            remapping(epoch,track_pair).median_replay_rate_diff=  remapping(epoch,track_pair).median_replay_rate_diff(index);
            % Difference median max instantaneous rate T1-T2 during REPLAY events where cell is active
            remapping(epoch,track_pair).mean_max_FR_replay_diff=  remapping(epoch,track_pair).mean_max_FR_replay_diff(index);
            
            % Difference mean FR T1-T2
  %          remapping(epoch,track_pair).replay_rate_diff=(remapping(epoch,track_pair).track1_mean_replay_rate-remapping(epoch,track_pair).track2_mean_replay_rate);
            remapping(epoch,track_pair).log_replay_rate_diff=(log(remapping(epoch,track_pair).track1_median_replay_rate)-log(remapping(epoch,track_pair).track2_median_replay_rate));
            % difference proportion
            remapping(epoch,track_pair).proportion_events_diff= (remapping(epoch,track_pair).proportion_events_cell_active_track1 - remapping(epoch,track_pair).proportion_events_cell_active_track2)';
            
            remapping(epoch,track_pair).tracks_compared = tracks_to_compare;
            
            % keep spikes for distribution
            remapping(epoch,track_pair).spikes1= remapping(epoch,track_pair).spikes1(index,:);
            remapping(epoch,track_pair).spikes2= remapping(epoch,track_pair).spikes2(index,:);
            remapping(epoch,track_pair).rate1= remapping(epoch,track_pair).rate1(index,:);
            remapping(epoch,track_pair).rate2= remapping(epoch,track_pair).rate2(index,:);
            
        end
    end

end

% switch method
%     case 'wcorr'
%         save rate_remapping_analysis_TRACK_PAIRS_wcorr remapping remapping_raw
%     case 'spearman'
%         save rate_remapping_analysis_TRACK_PAIRS_spearman remapping remapping_raw
%     case 'control_fixed_spike'
%         save rate_remapping_analysis_TRACK_PAIRS_wcorr_FIXED remapping remapping_raw
%     otherwise
%         save rate_remapping_analysis_TRACK_PAIRS remapping remapping_raw
% end

end