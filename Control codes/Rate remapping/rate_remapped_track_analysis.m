function  rate_remapped_track_analysis(current_folder,RExp,varargin)
% varargin: 
%           varargin{1}: seed, 
%           varargin{2}: 'shuffle_rates' or 'replay_analysis' to only
%           redo parts, otherwise re-runs everything

if exist('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\CONTROLS\rate_remapped')~=7
    mkdir('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\CONTROLS\rate_remapped');
end

if ~contains(current_folder,'\') && exist(['..\CONTROLS\rate_remapped\' current_folder])~=7 %when current folder is just the name of the folder
    mkdir(['..\CONTROLS\rate_remapped\' current_folder]);
elseif exist(current_folder)~=7 %for when current folder is the whole path
    mkdir(current_folder);
end


if ~isempty(varargin) && length(varargin)>1
    if strcmp(varargin{2},'shuffle_rates')
        if ~contains(current_folder,'\')
            %copy relevant data from main folder to chimeric folder
            copyfile('extracted_position.mat',['..\CONTROLS\rate_remapped\' current_folder]);
            copyfile('extracted_replay_events.mat',['..\CONTROLS\rate_remapped\' current_folder]);
            copyfile('extracted_sleep_state.mat',['..\CONTROLS\rate_remapped\' current_folder]);
            copyfile('extracted_clusters.mat',['..\CONTROLS\rate_remapped\' current_folder]);
            copyfile('extracted_place_fields_BAYESIAN.mat',['..\CONTROLS\rate_remapped\' current_folder]);
            copyfile('extracted_place_fields.mat',['..\CONTROLS\rate_remapped\' current_folder]);
            copyfile('extracted_waveforms.mat',['..\CONTROLS\rate_remapped\' current_folder]);
            
            cd(['..\CONTROLS\rate_remapped\' current_folder]);
        else
            %copy relevant data from main folder to chimeric folder
            idx = strfind(current_folder,'\');
            fname = current_folder;
            fname = current_folder(idx(end)+1:end);
            copyfile(['..\..\..\' fname '\extracted_position.mat'],current_folder);
            copyfile(['..\..\..\' fname '\extracted_replay_events.mat'],current_folder);
            copyfile(['..\..\..\' fname '\extracted_sleep_state.mat'],current_folder);
            copyfile(['..\..\..\' fname '\extracted_clusters.mat'],current_folder);
            copyfile(['..\..\..\' fname '\extracted_place_fields_BAYESIAN.mat'],current_folder);
            copyfile(['..\..\..\' fname '\extracted_place_fields.mat'],current_folder);
            copyfile(['..\..\..\' fname '\extracted_waveforms.mat'],current_folder);
            
            cd(current_folder)
        end

        create_rate_remapped_track('BAYESIAN',varargin{1});
        create_rate_remapped_track('NORMAL',varargin{1});
        remapping_pipeline('REPLAY',RExp);  %run batch_analysis for decoding replay and running shuffles

    elseif strcmp(varargin{2},'replay_analysis')
        if ~contains(current_folder,'\') 
            cd(['..\CONTROLS\rate_remapped\' current_folder]);
        else
            cd(current_folder)
        end
        remapping_pipeline('REPLAY',RExp);  %run batch_analysis for decoding replay and running shuffles
        %remapping_pipeline('ANALYSIS',RExp);
    end
        
else 
        %copy relevant data from main folder to chimeric folder
        copyfile('extracted_position.mat',['..\CONTROLS\rate_remapped\' current_folder]);
        copyfile('extracted_replay_events.mat',['..\CONTROLS\rate_remapped\' current_folder]);
        copyfile('extracted_sleep_state.mat',['..\CONTROLS\rate_remapped\' current_folder]);
        copyfile('extracted_clusters.mat',['..\CONTROLS\rate_remapped\' current_folder]);
        copyfile('extracted_place_fields_BAYESIAN.mat',['..\CONTROLS\rate_remapped\' current_folder]);
        copyfile('extracted_place_fields.mat',['..\CONTROLS\rate_remapped\' current_folder]);
        copyfile('extracted_waveforms.mat',['..\CONTROLS\rate_remapped\' current_folder]);
        
        if ~contains(current_folder,'\') 
            cd(['..\CONTROLS\rate_remapped\' current_folder]);
        else
            cd(current_folder)
        end 

        create_rate_remapped_track('BAYESIAN',varargin);
        create_rate_remapped_track('NORMAL',varargin);
        
        remapping_pipeline('REPLAY',RExp);  %run batch_analysis for decoding replay and running shuffles
        remapping_pipeline('ANALYSIS',RExp);
       
end

end