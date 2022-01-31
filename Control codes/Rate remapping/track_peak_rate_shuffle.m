% TRACK PEAK RATE SHUFFLE 
% Scales peak FR of place cells. Runs two types of controls based on rate change:
    % Rate intrinsic bias: runs control for intrinsic bias & formerly rate
        % remapping shuffle. The former uses real events but plots shuffled
        % peak FR in the track. The latter plots all shuffled data.
    % Rate detection control: After running the rate intrinsic bias, saves the indices of the "new detected events" from 
        % the shuffled data and uses them to plot the real events
% INPUT:
    % 'rate_intrinsic_bias', 'rate_detection_control', or 'rate_intrinsic_bias_distribution'. 'Rate_intrinsic_bias' needs
    % to have been run for being able to run 'rate_detection_control'.

function track_peak_rate_shuffle(folders,varargin)

% Set parameters and paths
save_path = [pwd '\CONTROLS\rate_remapped'];
if ~exist(save_path)
    mkdir(save_path);
end

parameters = list_of_parameters;
num_shuffles = 1000;
if isempty(folders)
    if exist('folders_to_process_remapping.mat')
        load([pwd '\folders_to_process_remapping.mat']);
    elseif exist('folders_to_process_REPLAY.mat')
        load([pwd '\folders_to_process_REPLAY.mat']);
    end
end
if contains(pwd,'Reward')
    folders= arrayfun(@(x) [pwd '\' folders{x}],1:length(folders),'UniformOutput',0);
    fd= strsplit(folders{1},'\');
    if length(fd)>1
        ctrl_folders= arrayfun(@(x) strsplit(folders{x},'\'),1:length(folders),'UniformOutput',0);
        ctrl_folders= arrayfun(@(x) ctrl_folders{x}{end},1:length(folders),'UniformOutput',0)';
    end
    shuffle_folders = arrayfun(@(x) [save_path '\' ctrl_folders{x}],1:length(folders),'UniformOutput',0);
else
    shuffle_folders= folders(:,1);
    ctrl_folders= folders(:,1);
end

if ~isempty(varargin) && strcmp(varargin{1},'rate_intrinsic_bias')
    %%% Intrinsic bias control - negative control
    % Use original events vs shuffled peak FR in track
%     if contains(pwd,'remapping')
%         RExp = {[],1,[],[],[]};
%     elseif contains(pwd,'Reward')
%         RExp= {[],[],[],[],[],[],...
%             [],[],[],[],[],[],...
%             [],[],[],[],[],[]};
%     end
    RExp= folders(:,2);
    for ses = 1 : length(shuffle_folders)
        cd([folders{ses}])
        disp(ctrl_folders{ses})
        rate_remapped_track_analysis(ctrl_folders{ses},RExp{ses},[],'shuffle_rates');
        sort_replay_events([],'wcorr');
        cd(save_path);
    end
    
    rate_remapping_TRACK_PAIRS(ctrl_folders,'wcorr',1);    %Runs former rate remapping 
    rate_remapping_TRACK_PAIRS(ctrl_folders,'rate_intrinsic_bias_control',1);    % Runs intrinsic bias 

    cd ..\..
    % Former rate remapping
%     plot_rate_remapping_NEW('x_var',{'place_field_diff','track_mean_rate_diff','place_field_diff'},'y_var',...
%         {'replay_spike_diff_nonZero','median_replay_rate_diff','mean_max_FR_replay_diff'},'epochs',{'PRE','POST'},'control','rate_remapped')
     plot_rate_remapping_NEW('subset',{'stable cells laps'},'epochs',{'PRE','POST'},'control','rate_remapped')
    % Intrinsic bias
    plot_rate_remapping_NEW('x_var',{'place_field_diff','track_mean_rate_diff','place_field_diff'},'y_var',...
        {'replay_spike_diff_nonZero','median_replay_rate_diff','mean_max_FR_replay_diff'},'epochs',{'PRE','POST'},...
        'control','rate_intrinsic_bias_control')
    
elseif ~isempty(varargin) && strcmp(varargin{1},'rate_detection_control')
    %%% Detection control - positive control
    % Uses indices from sign events detected with the shuffled data to plot real events vs real peak FR in track
    
    for ses = 1 : length(shuffle_folders)
        cd(shuffle_folders{ses})
        disp(shuffle_folders{ses})
        number_of_significant_replays(.05,3, 'wcorr', [],'rate_detection_control');
        sort_replay_events([],'rate_detection_control');
%         compare_replay_detection([],'rate_detection_control')
    end
    
    cd(save_path)
    rate_remapping_TRACK_PAIRS(shuffle_folders,'rate_detection_control',1);
    cd ..\..
    plot_rate_remapping_NEW('x_var',{'place_field_diff'},'y_var',{'mean_max_FR_replay_diff'},'y_label','Peak replay rate diff (Hz)',...
    'x_label','Peak in-field rate diff(Hz)','epochs',{'PRE','POST'},'control','rate_detection_control','subset','stable cells laps')


elseif ~isempty(varargin) && strcmp(varargin{1},'rate_intrinsic_bias_distribution')
    
    % Create shuffle distribution of Rs for 'Intrinsic Bias control'
    temp_folders = arrayfun(@(x) [save_path '\TEMP_' ctrl_folders{x}],1:length(folders),'UniformOutput',0);
    
    for s = 1 : num_shuffles
        s
        for ses = 1 : length(temp_folders)
            if ~exist(temp_folders{ses})
                mkdir(temp_folders{ses})
                cd(temp_folders{ses})
                copyfile([folders{ses} '\extracted_clusters.mat'],temp_folders{ses});
                copyfile([folders{ses} '\extracted_position.mat'],temp_folders{ses});
                copyfile([folders{ses} '\extracted_waveforms.mat'],temp_folders{ses});
                copyfile([folders{ses} '\extracted_place_fields_BAYESIAN.mat'],temp_folders{ses});
        
%                 copyfile(['..\..\..\' folders{ses} '\extracted_clusters.mat'],temp_folders{ses});
%                 copyfile(['..\..\..\' folders{ses} '\extracted_position.mat'],temp_folders{ses});
%                 copyfile(['..\..\..\' folders{ses} '\extracted_waveforms.mat'],temp_folders{ses});
%                 copyfile(['..\..\..\' folders{ses} '\extracted_place_fields_BAYESIAN.mat'],temp_folders{ses});
            else
                cd(temp_folders{ses})
                  copyfile([folders{ses} '\extracted_place_fields_BAYESIAN.mat'],temp_folders{ses});
%                 copyfile(['..\..\..\' folders{ses} '\extracted_place_fields_BAYESIAN.mat'],temp_folders{ses});
            end
            create_rate_remapped_track('BAYESIAN');%,parameters.rng_seed_remapping(ses)
        end
        cd(save_path)
        % Runs intrinsic bias
        [remapping, ~] = rate_remapping_TRACK_PAIRS(temp_folders,'rate_intrinsic_bias_control',0);
        [Fstat(s,:),pval(s,:)] = get_corr_value(remapping,'x_var',{'place_field_diff'},'y_var',....
        {'mean_max_FR_replay_diff'},'epochs',{'PRE','POST'});
    end
    
    % save variable
    intrinsic_rate_bias_shuffle_dist.Fstat = Fstat; %first column is PRE, second is POST
    intrinsic_rate_bias_shuffle_dist.pval = pval;
    save([save_path '\intrinsic_rate_bias_shuffle_dist.mat'],'intrinsic_rate_bias_shuffle_dist','-v7.3')
    
    % Delete temporary folders
    arrayfun(@(x) rmdir(temp_folders{x},'s'),1:length(temp_folders))
    
end

end


function [f,pval] = get_corr_value(remapping,varargin)

load([pwd '\Tables\subsets_of_cells.mat']);
available_scoring= {'wcorr','spearman'};
available_subsets= unique(subset_of_cells.subset)';

p= inputParser;
addParameter(p,'scoring','wcorr',@(x) ismember(x,available_scoring));
addParameter(p,'epochs',{'PRE','POST'},@iscell);
addParameter(p,'x_var',{'place_field_diff'},@(x) ischar(x) || iscell(x));
addParameter(p,'y_var',{'mean_max_FR_replay_diff'},@(x) ischar(x) || iscell(x));
addParameter(p,'common_cells_toggle',0,@ismatrix);
addParameter(p,'subset',[],@(x) ismember(x,available_subsets));
parse(p,varargin{:});

if ~isempty(p.Results.subset)
    subset= subset_of_cells.cell_IDs{strcmp(subset_of_cells.subset,p.Results.subset)};
else
    subset= [];
end
epochs= p.Results.epochs;
[~,epoch_idx]= ismember(epochs,vertcat(remapping.epoch));
x_var= p.Results.x_var;
y_var= p.Results.y_var;
if ischar(x_var)
    x_var= {x_var};
end
if ischar(y_var)
    x_var= {y_var};
end
parameters=list_of_parameters;

for track_pair = 1:size(remapping,2)
    
    % make sure repetitions are not discarded
    common_epoch_cells= multintersect(remapping(epoch_idx,track_pair).new_ID);
    
    for this_epoch=1:length(epochs)
        for this_var_pair=1:length(x_var)
            % check the var names are fields from loaded filename
            if ~isfield(remapping,x_var{this_var_pair}) ||  ~isfield(remapping,y_var{this_var_pair})
                error('wrong variable name entered, no plot created');
            end
            
            % select cells
            if p.Results.common_cells_toggle % keep only cells common to selected epochs
                epoch_cells = intersect(remapping(epoch_idx(this_epoch),track_pair).new_ID,common_epoch_cells);
            else % take all cells
                epoch_cells= remapping(epoch_idx(this_epoch),track_pair).new_ID;
            end
            
            % Subsets
            if ~isempty(subset)
                epoch_cells= intersect(epoch_cells,subset);
                common_subset_cells= multintersect(remapping(epoch_idx(this_epoch),track_pair).new_ID,epoch_cells,subset);
            else % if no subset, just put cells common to epochs as default
                common_subset_cells = intersect(remapping(epoch_idx(this_epoch),track_pair).new_ID,epoch_cells);
            end
            
            % maybe not most straightforward way but makes sense
            epoch_cells_ind= find(ismember(epoch_cells,remapping(epoch_idx(this_epoch),track_pair).new_ID));
            common_subset_cells_ind= find(ismember(common_subset_cells,remapping(epoch_idx(this_epoch),track_pair).new_ID));
            
            % Get correlation and p-value
            index_non_NaNs= intersect(find(~isnan(remapping(epoch_idx(this_epoch),track_pair).(x_var{this_var_pair}))),common_subset_cells_ind);
            if ~isempty(index_non_NaNs)
                lm = fitlm(remapping(epoch_idx(this_epoch),track_pair).(x_var{this_var_pair})(index_non_NaNs),...
                    remapping(epoch_idx(this_epoch),track_pair).(y_var{this_var_pair})(index_non_NaNs),'linear');
                [pval(this_epoch,this_var_pair),f(this_epoch,this_var_pair),~] = coefTest(lm);

            end
        end
    end
end
end
