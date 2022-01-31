function replay_rate_shuffle_analysis(folders,method)
% works in separate folder (CONTROLS\fixed_spike_count)
% input list of folders to process
% and vector of reexposure options

master_folder = pwd;
shuffle_folder= 'X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\CONTROLS\replay_rate_shuffle';
re_exposure= folders(:,2);
folders= folders(:,1);
parameters= list_of_parameters;

    switch method
        case 'replay_rate_shuffle'
            for this_session = 1 : size(folders,1)

                    cd([master_folder '\' folders{this_session}])
                    disp(folders{this_session})
                    RExp = re_exposure{this_session};
                    current_folder= folders{this_session};

                    if exist(['..\CONTROLS\replay_rate_shuffle\' current_folder])~=7
                        mkdir(['..\CONTROLS\replay_rate_shuffle\' current_folder]);
                    end
                    % copy over files needed 
                    copyfile('extracted_replay_events.mat',['..\CONTROLS\replay_rate_shuffle\' current_folder]);
                    copyfile('decoded_replay_events.mat',['..\CONTROLS\replay_rate_shuffle\' current_folder]);
                    copyfile('decoded_replay_events_segments.mat',['..\CONTROLS\replay_rate_shuffle\' current_folder]);
                    copyfile('extracted_clusters.mat',['..\CONTROLS\replay_rate_shuffle\' current_folder]);
                    copyfile('extracted_place_fields_BAYESIAN.mat',['..\CONTROLS\replay_rate_shuffle\' current_folder]);
                    copyfile('extracted_position.mat',['..\CONTROLS\replay_rate_shuffle\' current_folder]);
                    copyfile('replayEvents_bayesian_spike_count.mat',['..\CONTROLS\replay_rate_shuffle\' current_folder]);
                    copyfile('extracted_sleep_state.mat',['..\CONTROLS\replay_rate_shuffle\' current_folder]);
                    if exist('significant_replay_events_wcorr.mat')==2
                        copyfile('significant_replay_events_wcorr.mat',['..\CONTROLS\replay_rate_shuffle\' current_folder]);
                    elseif exist('significant_replay_events_wcorr_individual_exposures.mat')==2
                         copyfile('significant_replay_events_wcorr_individual_exposures.mat',['..\CONTROLS\replay_rate_shuffle\' current_folder]);
                    end

                    cd(['..\CONTROLS\replay_rate_shuffle\' current_folder]);    

                    replay_decoding('replay_rate_shuffle');
                    % run shuffles
                    remapping_pipeline('SHUFFLES',RExp);

                    %get significant replay events
                    if isempty(RExp) % no reexposure
                        number_of_significant_replays(0.05,3,'wcorr',[]); % pval, ripple zscore, method, reexposure
                        sort_replay_events([],'wcorr'); 
                    else
                        % reexposure
                        number_of_significant_replays(0.05,3,'wcorr',2);
                        sort_replay_events([1 3;2 4],'wcorr'); 
                    end
            end
            cd(shuffle_folder);
            rate_remapping_TRACK_PAIRS(folders,'replay_rate_shuffle_control');
            rate_remapping_TRACK_PAIRS(folders,'replay_rate_shuffle_detection');
            cd(master_folder);
            plot_rate_remapping_NEW('epochs',{'PRE','POST'},'control','replay_rate_shuffle_control','subset','stable cells laps',...
                'x_label','Peak Rate Change (Hz)','y_label','Replay Rate Change (Hz)')
            plot_rate_remapping_NEW('epochs',{'PRE','POST'},'control','replay_rate_shuffle_detection','subset','stable cells laps',...
                'x_label','Peak Rate Change (Hz)','y_label','Replay Rate Change (Hz)')
            
        case 'replay_rate_shuffles_only'
            % Set parameters and paths
            save_path = 'X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\CONTROLS\replay_rate_shuffle';
              % Create shuffle distribution of Rs - made consistent with
              % track_peak_rate_shuffle
                temp_folders = arrayfun(@(x) [shuffle_folder '\TEMP_' folders{x}],1:length(folders),'UniformOutput',0);
                Fstat= []; pval= [];
                for s = 1 : 1000
                    s
                    tic;
                    for this_session = 1 : length(temp_folders)
                        if ~exist(temp_folders{this_session})
                            mkdir(temp_folders{this_session})
                        end
                        cd(temp_folders{this_session});
%                         disp(temp_folders{this_session})
                        copyfile(['..\..\..\' folders{this_session} '\replayEvents_bayesian_spike_count.mat'],temp_folders{this_session});
                        copyfile(['..\..\..\' folders{this_session} '\extracted_place_fields_BAYESIAN.mat'],temp_folders{this_session});
                        copyfile(['..\..\..\' folders{this_session} '\decoded_replay_events.mat'],temp_folders{this_session});
                        if exist(['..\..\..\' folders{this_session} '\significant_replay_events_wcorr.mat'])==2
                            copyfile(['..\..\..\' folders{this_session} '\significant_replay_events_wcorr.mat'],temp_folders{this_session});
                        elseif exist(['..\..\..\' folders{this_session} '\significant_replay_events_wcorr_individual_exposures.mat'])==2
                             copyfile(['..\..\..\' folders{this_session} '\significant_replay_events_wcorr_individual_exposures.mat'],temp_folders{this_session});
                        end
                    end 
                    parfor this_session = 1 : length(temp_folders)
                        cd(temp_folders{this_session});
                        shuffle_spike_rate_replay;
                    end
                    
                    cd(shuffle_folder)
                    % Runs intrinsic bias
                    [remapping, ~] = rate_remapping_TRACK_PAIRS(temp_folders,'replay_rate_shuffle_control',0);
                    cd(master_folder)
                    [pval(s,:),Fstat(s,:)] = plot_rate_remapping_NEW('use_mat',remapping,'x_var',{'place_field_diff'},'y_var',....
                    {'mean_max_FR_replay_diff'},'epochs',{'PRE','POST'});
                    close all;
                    
                    % save every 50 iterations
                    if mod(s,50)==0
                        intrinsic_replay_rate_bias_shuffle_dist.Fstat = Fstat; %first column is PRE, second is POST
                        intrinsic_replay_rate_bias_shuffle_dist.pval = pval;
                        save([save_path '\intrinsic_replay_rate_bias_shuffle_dist.mat'],'intrinsic_replay_rate_bias_shuffle_dist','-v7.3');
                        disp('saved file')
                    end
                    disp('running one shuffle took..')
                    toc
                end

                    % save variable
                    intrinsic_replay_rate_bias_shuffle_dist.Fstat = Fstat; %first column is PRE, second is POST
                    intrinsic_replay_rate_bias_shuffle_dist.pval = pval;
                    save([save_path '\intrinsic_replay_rate_bias_shuffle_dist.mat'],'intrinsic_replay_rate_bias_shuffle_dist','-v7.3')

                    % Delete temporary folders
                    arrayfun(@(x) rmdir(temp_folders{x},'s'),1:length(temp_folders))
            
    end




end