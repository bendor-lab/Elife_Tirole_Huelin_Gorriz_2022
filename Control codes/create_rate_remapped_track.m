function create_rate_remapped_track(option)

parameters = list_of_parameters;

load('extracted_clusters.mat');
load('extracted_position.mat');
if exist('extracted_waveforms.mat','file')
    load('extracted_waveforms.mat');
else
    disp('no extracted_waveforms.mat file');
    allclusters_waveform = [];
end

current_directory = pwd;
parameters = list_of_parameters;

if strcmp(option,'BAYESIAN')
    load extracted_place_fields_BAYESIAN;
    place_fields = place_fields_BAYESIAN;
    clear place_fields_BAYESIAN;
    w = [1 1];
    x_bins_width=parameters.x_bins_width_bayesian;
else
    load extracted_place_fields;
    w = gausswin(parameters.place_field_smoothing);
    x_bins_width=parameters.x_bins_width;
end
%for smoothing
w = w./sum(w);

number_of_rate_remapped_tracks = length(place_fields.track);
number_of_tracks = length(place_fields.track);
number_of_neurons = length(place_fields.track(1).peak);

% Run threshold on pyramidal cells: half-width amplitude
if ~isempty(allclusters_waveform)
    PC_indices = [allclusters_waveform.half_width] > parameters.half_width_threshold; % cells that pass treshold of pyramidal cell half width
    pyramidal_cells = [allclusters_waveform(PC_indices).converted_ID];
end

for track_id = 1 : number_of_rate_remapped_tracks
    
    place_fields_RATE_REMAPPED.track(track_id).time_window = [NaN NaN];
    place_fields_RATE_REMAPPED.track(track_id).x_bin_centres = place_fields.track(track_id).x_bin_centres;
    place_fields_RATE_REMAPPED.track(track_id).x_bin_edges = place_fields.track(track_id).x_bin_edges;
    place_fields_RATE_REMAPPED.track(track_id).x_bins_width = place_fields.track(track_id).x_bins_width;
    place_fields_RATE_REMAPPED.track(track_id).dwell_map = [];
    place_fields_RATE_REMAPPED.track(track_id).good_cells = [];
    place_fields_RATE_REMAPPED.track(track_id).good_cells_LIBERAL =[];
    place_fields_RATE_REMAPPED.track(track_id).sorted =[];
    place_fields_RATE_REMAPPED.track(track_id).sorted_good_cells =[];
    place_fields_RATE_REMAPPED.track(track_id).sorted_good_cells_LIBERAL =[];
    place_fields_RATE_REMAPPED.track(track_id).unique_cells =[];
    
    % Sort good cells by raw peak firing rate
    [raw_peak_distribution,I] = sort(place_fields.track(track_id).raw_peak(place_fields.track(track_id).good_cells));
    fraction_of_cells = (0:(length(place_fields.track(track_id).good_cells)-1))/(length(place_fields.track(track_id).good_cells)-1);
    
    for j=1:number_of_neurons
        
        random_peak_rate = interp1(fraction_of_cells,raw_peak_distribution,rand_FR,'linear');
        
        if isempty(find(place_fields.track(track_id).good_cells==j))
            scaling_factor=1;
        else
            scaling_factor=random_peak_rate./place_fields.track(track_id).raw_peak(j);
        end
        
        place_fields_RATE_REMAPPED.track(track_id).spike_hist{j}=scaling_factor*place_fields.track(track_id).spike_hist{j};
        place_fields_RATE_REMAPPED.track(track_id).raw{j}=scaling_factor*place_fields.track(track_id).raw{j};
        
        % Get place field information
        place_fields_RATE_REMAPPED.track(track_id).smooth{j}         =  filtfilt(w,1,place_fields_RATE_REMAPPED.track(track_id).raw{j}); %smooth pl field
        place_fields_RATE_REMAPPED.track(track_id).centre_of_mass(j) =  sum(place_fields_RATE_REMAPPED.track(track_id).smooth{j}.* place_fields_RATE_REMAPPED.track(track_id).x_bin_centres)/sum(place_fields_RATE_REMAPPED.track(track_id).smooth{j});  %averaged center
        [place_fields_RATE_REMAPPED.track(track_id).peak(j) , index] =  max(place_fields_RATE_REMAPPED.track(track_id).smooth{j}); %peak of smoothed place field and index of peak (center)
        if place_fields_RATE_REMAPPED.track(track_id).peak(j) ~=0
            if length(index)>1 % very rare exception where you have multiple peaks of same height....
                index= index(1);
            end
            place_fields_RATE_REMAPPED.track(track_id).centre(j) =  place_fields_RATE_REMAPPED.track(track_id).x_bin_centres(index);
            
        else
            place_fields_RATE_REMAPPED.track(track_id).centre(j) = NaN;
        end
        place_fields_RATE_REMAPPED.track(track_id).raw_peak(j)          = max(place_fields_RATE_REMAPPED.track(track_id).raw{j}); % raw pl field peak
        place_fields_RATE_REMAPPED.track(track_id).mean_rate_session(j) = NaN; %mean firing rate
        place_fields_RATE_REMAPPED.track(track_id).mean_rate_track(j)   = sum(place_fields_RATE_REMAPPED.track(track_id).spike_hist{j})/(place_fields_RATE_REMAPPED.track(track_id).time_window(2)-place_fields_RATE_REMAPPED.track(track_id).time_window(1));
        if place_fields_RATE_REMAPPED.track(track_id).peak(j) ~=0
            place_fields_RATE_REMAPPED.track(track_id).half_max_width(j) = x_bins_width*half_max_width(place_fields_RATE_REMAPPED.track(track_id).smooth{j}); %finds half width of smoothed place field (width from y values closest to 50% of peak)
        else
            place_fields_RATE_REMAPPED.track(track_id).half_max_width(j) = NaN;
        end
        
    end
    
    %calculate skagges information
    place_fields_RATE_REMAPPED.track(track_id).skaggs_info= place_fields.track(track_id).skaggs_info;
    place_fields_RATE_REMAPPED.track(track_id).good_cells = place_fields.track(track_id).good_cells; % Check that the cells that passed the threshold are pyramidal cells
    place_fields_RATE_REMAPPED.track(track_id).good_cells_LIBERAL =  place_fields.track(track_id).good_cells_LIBERAL; % Check that the cells that passed the threshold are pyramidal cells
    place_fields_RATE_REMAPPED.track(track_id).sorted =  place_fields.track(track_id).sorted;
    place_fields_RATE_REMAPPED.track(track_id).sorted_good_cells = place_fields.track(track_id).sorted_good_cells;
    place_fields_RATE_REMAPPED.track(track_id).sorted_good_cells_LIBERAL = place_fields.track(track_id).sorted_good_cells_LIBERAL;
end

%% Classify cells as good place cells, interneuron, pyramidal cells & other cells

place_fields_RATE_REMAPPED.good_place_cells = place_fields.good_place_cells;
place_fields_RATE_REMAPPED.good_place_cells_LIBERAL = place_fields.good_place_cells_LIBERAL;
place_fields_RATE_REMAPPED.unique_cells = place_fields.unique_cells;  % all cells that have good place fields only on a single track
place_fields_RATE_REMAPPED.interneurons = place_fields.interneurons;
place_fields_RATE_REMAPPED.pyramidal_cells = place_fields.pyramidal_cells;
place_fields_RATE_REMAPPED.all_cells= place_fields.all_cells;
place_fields_RATE_REMAPPED.other_cells = place_fields.other_cells;


if exist('rate_remapped')~=7
    mkdir rate_remapped;
end
cd rate_remapped;

if strcmp(option,'BAYESIAN')
    place_fields_BAYESIAN=place_fields_RATE_REMAPPED;
    save extracted_place_fields_BAYESIAN place_fields_BAYESIAN;
else
    place_fields=place_fields_RATE_REMAPPED;
    save extracted_place_fields place_fields;
end

cd(current_directory);

end




function half_width = half_max_width(place_field)
%interpolate place field to get better resolution
new_step_size=0.1;  %decrease value to get finer resolution interpolation of place field
place_field_resampled=interp1(1:length(place_field),place_field,1:new_step_size:length(place_field),'linear');
[peak,index] = max(place_field_resampled); %finds smoothed place field peak firing rate (FR)
for i = index : length(place_field_resampled)
    if place_field_resampled(i)<peak/2 %finds the point after the peak where the FR is half the peak FR
        break;
    end
end
for j = index : -1 : 1 %finds the point before the peak where the FR is half the peak FR
    if place_field_resampled(j)<peak/2
        break;
    end
end
half_width = new_step_size*(i-j); %distance between half-peaks
%(calculated in indicies of original place field, but converted to distance in cm in function above)
end
