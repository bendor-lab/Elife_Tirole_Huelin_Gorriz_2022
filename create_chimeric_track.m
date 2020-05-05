function create_chimeric_track(option)
current_directory=pwd;
parameters=list_of_parameters;
if strcmp(option,'BAYESIAN')
    load extracted_place_fields_BAYESIAN;
    place_fields=place_fields_BAYESIAN;
    clear place_fields_BAYESIAN;
    w= [1 1];
else
    load extracted_place_fields;
     w= gausswin(parameters.place_field_smoothing);
end
w=w./sum(w);  %for smoothing

number_of_chimeric_tracks=length(place_fields.track);
number_of_tracks=length(place_fields.track);
number_of_neurons=length(place_fields.track(1).peak);

 
for k=1:number_of_chimeric_tracks
    place_fields_CHIMERIC.track(k).time_window=[NaN NaN];
    place_fields_CHIMERIC.track(k).x_bin_centres=place_fields.track(1).x_bin_centres;
    place_fields_CHIMERIC.track(k).x_bin_edges=place_fields.track(1).x_bin_edges;
    place_fields_CHIMERIC.track(k).x_bins_width=place_fields.track(1).x_bins_width;
    place_fields_CHIMERIC.track(k).dwell_map=[];
    place_fields_CHIMERIC.track(k).good_cells=[];
    place_fields_CHIMERIC.track(k).good_cells_LIBERAL=[];
    place_fields_CHIMERIC.track(k).sorted=[];
    place_fields_CHIMERIC.track(k).sorted_good_cells=[];
    place_fields_CHIMERIC.track(k).sorted_good_cells_LIBERAL=[];
    place_fields_CHIMERIC.track(k).unique_cells=[];
    
    for j=1:number_of_neurons
        random_track=randperm(number_of_tracks);
        random_track=random_track(1);
        reflection=randperm(2); reflection=reflection(1);
        if reflection==1
            place_fields_CHIMERIC.track(k).spike_hist{j}=place_fields.track(random_track).spike_hist{j};
            place_fields_CHIMERIC.track(k).raw{j}=place_fields.track(random_track).raw{j};
        elseif reflection==2
            place_fields_CHIMERIC.track(k).spike_hist{j}=fliplr(place_fields.track(random_track).spike_hist{j});
            place_fields_CHIMERIC.track(k).raw{j}=fliplr(place_fields.track(random_track).raw{j});
        end
        
        
        place_fields_CHIMERIC.track(k).smooth{j}         = filtfilt(w,1,place_fields_CHIMERIC.track(k).raw{j}); %smooth pl field
        place_fields_CHIMERIC.track(k).centre_of_mass(j) = sum(place_fields_CHIMERIC.track(k).smooth{j}.*place_fields_CHIMERIC.track(k).x_bin_centres/sum(place_fields_CHIMERIC.track(k).smooth{j}));  %averaged center
        [place_fields_CHIMERIC.track(k).peak(j) , index] = max(place_fields_CHIMERIC.track(k).smooth{j}); %peak of smoothed place field and index of peak (center)
        place_fields_CHIMERIC.track(k).centre(j) =  place_fields_CHIMERIC.track(k).x_bin_centres(index);
        place_fields_CHIMERIC.track(k).raw_peak(j)          = max(place_fields_CHIMERIC.track(k).raw{j}); % raw pl field peak
        
        place_fields_CHIMERIC.track(k).mean_rate_session(j)=place_fields.track(random_track).mean_rate_session(j);
        place_fields_CHIMERIC.track(k).mean_rate_track(j)=place_fields.track(random_track).mean_rate_track(j);
        place_fields_CHIMERIC.track(k).half_max_width(j)=place_fields.track(random_track).half_max_width(j);

        place_fields_CHIMERIC.track(k).skaggs_info(j)=place_fields.track(random_track).skaggs_info(j);
        place_fields_CHIMERIC.track(k).non_visited_bins=place_fields.track(random_track).non_visited_bins;

        if length(find(place_fields.track(random_track).good_cells==j))==1
            place_fields_CHIMERIC.track(k).good_cells=[place_fields_CHIMERIC.track(k).good_cells j];
        end
        if length(find(place_fields.track(random_track).good_cells_LIBERAL==j))==1
            place_fields_CHIMERIC.track(k).good_cells_LIBERAL=[place_fields_CHIMERIC.track(k).good_cells_LIBERAL j];
        end
    end

    [~,index] = sort(place_fields_CHIMERIC.track(k).centre);
    place_fields_CHIMERIC.track(k).sorted = index;
    [~,index] = sort(place_fields_CHIMERIC.track(k).centre(place_fields_CHIMERIC.track(k).good_cells));
    place_fields_CHIMERIC.track(k).sorted_good_cells = place_fields_CHIMERIC.track(k).good_cells(index);
    [~,index] = sort(place_fields_CHIMERIC.track(k).centre(place_fields_CHIMERIC.track(k).good_cells_LIBERAL));
    place_fields_CHIMERIC.track(k).sorted_good_cells_LIBERAL = place_fields_CHIMERIC.track(k).good_cells_LIBERAL(index);
end




%% Classify cells as good place cells, interneuron, pyramidal cells & other cells
good_place_cells=[]; track=[];
for track_id=1:number_of_tracks %good cells classfication
    good_place_cells = [good_place_cells place_fields_CHIMERIC.track(track_id).sorted_good_cells];
    track =[track track_id*ones(size(place_fields_CHIMERIC.track(track_id).sorted_good_cells))];
end
place_fields_CHIMERIC.good_place_cells = unique(good_place_cells);

good_place_cells_LIBERAL=[];
for track_id=1:number_of_tracks %good cells classfication
    good_place_cells_LIBERAL = [good_place_cells_LIBERAL place_fields_CHIMERIC.track(track_id).sorted_good_cells_LIBERAL];
end
place_fields_CHIMERIC.good_place_cells_LIBERAL = unique(good_place_cells_LIBERAL);

% cells that are unique for each track
unique_cells=[]; 
for track_id = 1:number_of_tracks
    place_fields_CHIMERIC.track(track_id).unique_cells = setdiff(good_place_cells(track==track_id),good_place_cells(track~=track_id),'stable');
    unique_cells = [unique_cells, place_fields_CHIMERIC.track(track_id).unique_cells];
end
place_fields_CHIMERIC.unique_cells = unique_cells;  % all cells that have good place fields only on a single track

%interneurons classfication
interneurons=[];
for track_id=1:number_of_tracks
    interneurons = [interneurons find(place_fields_CHIMERIC.track(track_id).mean_rate_track > parameters.max_mean_rate)];
end
place_fields_CHIMERIC.interneurons = unique(interneurons);

%putative pyramidal cells classification
putative_pyramidal_cells=[];
for track_id=1:number_of_tracks % putative pyramidal cells that pass the 'Pyramidal type' threshold (but not need to be place cells)
    putative_pyramidal_cells = [putative_pyramidal_cells find(place_fields_CHIMERIC.track(track_id).mean_rate_track <= parameters.max_mean_rate)];
end


place_fields_CHIMERIC.pyramidal_cells = place_fields.pyramidal_cells;
place_fields_CHIMERIC.all_cells=place_fields.all_cells;
other_cells = setdiff(place_fields.all_cells,good_place_cells,'stable'); %find the excluded putative pyramidal cells
place_fields_CHIMERIC.other_cells = setdiff(other_cells,interneurons,'stable'); %remove also the interneurons




%to compare real place fields on track one with chimeric track

place_fields_compare.track(1)=place_fields.track(1);
place_fields_compare.track(2)=place_fields_CHIMERIC.track(1);
plot_place_fields(place_fields_compare);



%to compare all chimeric tracks

%plot_place_fields(place_fields_CHIMERIC)

if exist('chimeric')~=7
    mkdir chimeric;
end
cd chimeric;

if strcmp(option,'BAYESIAN')
    place_fields_BAYESIAN=place_fields_CHIMERIC;
    save extracted_place_fields_BAYESIAN place_fields_BAYESIAN;
else
    place_fields=place_fields_CHIMERIC;
    save extracted_place_fields place_fields;
end

cd(current_directory);

end