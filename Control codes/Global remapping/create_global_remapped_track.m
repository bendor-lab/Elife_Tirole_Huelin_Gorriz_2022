% GLOBAL REMAPPING
% Shuffle good place fields within a track using cell IDs from good cells from any of the tracks.
% option is: 'bayesian' or [] in which case it will use normal place fields
% varargin: specify random seed or not
function create_global_remapped_track(option,varargin)

current_directory = pwd;
parameters = list_of_parameters;
if ~isempty(varargin)
    disp('using random seed specified..')
    rng(cell2mat(varargin{1}));
end
load extracted_clusters.mat
if strcmp(option,'BAYESIAN')
    load extracted_place_fields_BAYESIAN;
    place_fields = place_fields_BAYESIAN;
    clear place_fields_BAYESIAN;
else
    load extracted_place_fields;
end

number_of_chimeric_tracks = length(place_fields.track);
number_of_tracks = length(place_fields.track);

place_fields_GLOBAL_REMAPPED = place_fields;  %by default make a copy of all cells

for k = 1 : number_of_chimeric_tracks
    
    % Allocate variables
    place_fields_GLOBAL_REMAPPED.track(k).good_cells=[];
    place_fields_GLOBAL_REMAPPED.track(k).good_cells_LIBERAL=[];
    place_fields_GLOBAL_REMAPPED.track(k).sorted=[];
    place_fields_GLOBAL_REMAPPED.track(k).sorted_good_cells=[];
    place_fields_GLOBAL_REMAPPED.track(k).sorted_good_cells_LIBERAL=[];
    place_fields_GLOBAL_REMAPPED.track(k).unique_cells=[];
    
    % Random permutation of good place cells from all tracks
    random_cell_index = randperm(length(place_fields.good_place_cells));
    random_cell = place_fields.good_place_cells(random_cell_index);
    original_cell = place_fields.good_place_cells;
    
    place_fields_GLOBAL_REMAPPED.id_conversion(k).original_good_cells= original_cell;
    place_fields_GLOBAL_REMAPPED.id_conversion(k).permuted_good_cells= random_cell;
    if ~isempty(varargin)
      place_fields_GLOBAL_REMAPPED.seed= varargin{1};
    end
    
    % Create a new shuffled structure, where we swap good cells (each original cells is assigned to another random cell)
    % Swap can happen between tracks (e.g. good cell from track 1 is now in track 2)
    for j = 1 : length(random_cell) %only swap good cells
        place_fields_GLOBAL_REMAPPED.track(k).spike_hist{original_cell(j)} = place_fields.track(k).spike_hist{random_cell(j)};
        place_fields_GLOBAL_REMAPPED.track(k).raw{original_cell(j)} = place_fields.track(k).raw{random_cell(j)};
        
        place_fields_GLOBAL_REMAPPED.track(k).smooth{original_cell(j)}         = place_fields.track(k).smooth{random_cell(j)};
        place_fields_GLOBAL_REMAPPED.track(k).centre_of_mass(original_cell(j)) = place_fields.track(k).centre_of_mass(random_cell(j));
        place_fields_GLOBAL_REMAPPED.track(k).peak(original_cell(j)) =  place_fields.track(k).peak(random_cell(j));
        place_fields_GLOBAL_REMAPPED.track(k).centre(original_cell(j)) =  place_fields.track(k).centre(random_cell(j));
        place_fields_GLOBAL_REMAPPED.track(k).raw_peak(original_cell(j)) =  place_fields.track(k).raw_peak(random_cell(j));
        
        place_fields_GLOBAL_REMAPPED.track(k).mean_rate_session(original_cell(j)) = place_fields.track(k).mean_rate_session(random_cell(j));
        place_fields_GLOBAL_REMAPPED.track(k).mean_rate_track(original_cell(j)) = place_fields.track(k).mean_rate_track(random_cell(j));
        place_fields_GLOBAL_REMAPPED.track(k).half_max_width(original_cell(j))=place_fields.track(k).half_max_width(random_cell(j));
        
        place_fields_GLOBAL_REMAPPED.track(k).skaggs_info(original_cell(j)) = place_fields.track(k).skaggs_info(random_cell(j));
        place_fields_GLOBAL_REMAPPED.track(k).non_visited_bins = place_fields.track(k).non_visited_bins;
                
        % If the current random cell belongs to the good cells from this track,include in good cells for that track
        if length(find(place_fields.track(k).good_cells==random_cell(j)))==1
            place_fields_GLOBAL_REMAPPED.track(k).good_cells=[place_fields_GLOBAL_REMAPPED.track(k).good_cells original_cell(j)];
        end
        if length(find(place_fields.track(k).good_cells_LIBERAL==random_cell(j)))==1
            place_fields_GLOBAL_REMAPPED.track(k).good_cells_LIBERAL=[place_fields_GLOBAL_REMAPPED.track(k).good_cells_LIBERAL original_cell(j)];
        end
    end
    
    % Sort new cells by place field centre
    [~,index] = sort(place_fields_GLOBAL_REMAPPED.track(k).centre);
    place_fields_GLOBAL_REMAPPED.track(k).sorted = index;
    [~,index] = sort(place_fields_GLOBAL_REMAPPED.track(k).centre(place_fields_GLOBAL_REMAPPED.track(k).good_cells));
    place_fields_GLOBAL_REMAPPED.track(k).sorted_good_cells = place_fields_GLOBAL_REMAPPED.track(k).good_cells(index);
    [~,index] = sort(place_fields_GLOBAL_REMAPPED.track(k).centre(place_fields_GLOBAL_REMAPPED.track(k).good_cells_LIBERAL));
    place_fields_GLOBAL_REMAPPED.track(k).sorted_good_cells_LIBERAL = place_fields_GLOBAL_REMAPPED.track(k).good_cells_LIBERAL(index);
end




%% Classify cells as good place cells, interneuron, pyramidal cells & other cells
%interneurons classfication
interneurons = find(place_fields_GLOBAL_REMAPPED.mean_rate > parameters.max_mean_rate);
place_fields_GLOBAL_REMAPPED.interneurons=interneurons;

good_place_cells=[]; track=[];
for track_id=1:number_of_tracks %good cells classfication
    good_place_cells = [good_place_cells place_fields_GLOBAL_REMAPPED.track(track_id).sorted_good_cells];
    track =[track track_id*ones(size(place_fields_GLOBAL_REMAPPED.track(track_id).sorted_good_cells))];
end
place_fields_GLOBAL_REMAPPED.good_place_cells = unique(good_place_cells);

good_place_cells_LIBERAL=[];
for track_id=1:number_of_tracks %good cells classfication
    good_place_cells_LIBERAL = [good_place_cells_LIBERAL place_fields_GLOBAL_REMAPPED.track(track_id).sorted_good_cells_LIBERAL];
end
place_fields_GLOBAL_REMAPPED.good_place_cells_LIBERAL = unique(good_place_cells_LIBERAL);

% cells that are unique for each track
unique_cells=[]; 
for track_id = 1:number_of_tracks
    place_fields_GLOBAL_REMAPPED.track(track_id).unique_cells = setdiff(good_place_cells(track==track_id),good_place_cells(track~=track_id),'stable');
    unique_cells = [unique_cells, place_fields_GLOBAL_REMAPPED.track(track_id).unique_cells];
end
place_fields_GLOBAL_REMAPPED.unique_cells = unique_cells;  % all cells that have good place fields only on a single track


%putative pyramidal cells classification
putative_pyramidal_cells=[];
for track_id=1:number_of_tracks % putative pyramidal cells that pass the 'Pyramidal type' threshold (but not need to be place cells)
    putative_pyramidal_cells = [putative_pyramidal_cells find(place_fields_GLOBAL_REMAPPED.track(track_id).mean_rate_track <= parameters.max_mean_rate)];
end

place_fields_GLOBAL_REMAPPED.pyramidal_cells = place_fields.pyramidal_cells;
place_fields_GLOBAL_REMAPPED.all_cells =1:max(clusters.id_conversion(:,1));
other_cells = setdiff(1:max(clusters.id_conversion(:,1)),good_place_cells,'stable'); %find the excluded putative pyramidal cells
place_fields_GLOBAL_REMAPPED.other_cells = setdiff(other_cells,interneurons,'stable'); %remove also the interneurons

% SAFETY CHECK: to compare real place fields on track one with chimeric track
    % place_fields_compare.track(1)=place_fields.track(1);
    % place_fields_compare.track(2)=place_fields_GLOBAL_REMAPPED.track(1);
    % plot_place_fields(place_fields_compare);


%to compare all chimeric tracks

% if exist('global_remapped')~=7
%     mkdir global_remapped;
%     cd global_remapped;
% end

if strcmp(option,'BAYESIAN')
    place_fields_BAYESIAN=place_fields_GLOBAL_REMAPPED;
    save extracted_place_fields_BAYESIAN place_fields_BAYESIAN;
else
    place_fields=place_fields_GLOBAL_REMAPPED;
    save extracted_place_fields place_fields;
end

cd(current_directory);

end