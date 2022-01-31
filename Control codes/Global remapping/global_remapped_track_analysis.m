function  global_remapped_track_analysis(current_folder,RExp,varargin)
% varargin: 
%           varargin{1}: seed, 
%           varargin{2}: 'shuffle_track' or 'replay_analysis' to only
%           redo parts, otherwise re-runs everything

if exist('..\CONTROLS\global_remapped')~=7
    mkdir('..\CONTROLS\global_remapped');
end

if exist(['..\CONTROLS\global_remapped\' current_folder])~=7
    mkdir(['..\CONTROLS\global_remapped\' current_folder]);
end

if ~isempty(varargin) && length(varargin)>1
    if strcmp(varargin{2},'shuffle_track')
        %copy relevant data from main folder to chimeric folder
        copyfile('extracted_position.mat',['..\CONTROLS\global_remapped\' current_folder]);
        copyfile('extracted_replay_events.mat',['..\CONTROLS\global_remapped\' current_folder]);
        copyfile('extracted_sleep_state.mat',['..\CONTROLS\global_remapped\' current_folder]);
        copyfile('extracted_clusters.mat',['..\CONTROLS\global_remapped\' current_folder]);
        copyfile('extracted_place_fields_BAYESIAN.mat',['..\CONTROLS\global_remapped\' current_folder]);
        copyfile('extracted_place_fields.mat',['..\CONTROLS\global_remapped\' current_folder]);
        copyfile('extracted_waveforms.mat',['..\CONTROLS\global_remapped\' current_folder]);
        
        cd(['..\CONTROLS\global_remapped\' current_folder]); 

        create_global_remapped_track('BAYESIAN',varargin{1});
        create_global_remapped_track('NORMAL',varargin{1});
        
    elseif strcmp(varargin{2},'replay_analysis')
        cd(['..\CONTROLS\global_remapped\' current_folder]); 
        remapping_pipeline('REPLAY',RExp);  %run batch_analysis for decoding replay and running shuffles
        remapping_pipeline('ANALYSIS',RExp);
    end
        
else 
        %copy relevant data from main folder to chimeric folder
        copyfile('extracted_position.mat',['..\CONTROLS\global_remapped\' current_folder]);
        copyfile('extracted_replay_events.mat',['..\CONTROLS\global_remapped\' current_folder]);
        copyfile('extracted_sleep_state.mat',['..\CONTROLS\global_remapped\' current_folder]);
        copyfile('extracted_clusters.mat',['..\CONTROLS\global_remapped\' current_folder]);
        copyfile('extracted_place_fields_BAYESIAN.mat',['..\CONTROLS\global_remapped\' current_folder]);
        copyfile('extracted_place_fields.mat',['..\CONTROLS\global_remapped\' current_folder]);
        copyfile('extracted_waveforms.mat',['..\CONTROLS\global_remapped\' current_folder]);
        
        cd(['..\CONTROLS\global_remapped\' current_folder]); 

        create_global_remapped_track('BAYESIAN',varargin);
        create_global_remapped_track('NORMAL',varargin);
        
        remapping_pipeline('REPLAY',RExp);  %run batch_analysis for decoding replay and running shuffles
        remapping_pipeline('ANALYSIS',RExp);
       
end

end