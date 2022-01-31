function decoding_error= decoding_comparison(folders,varargin)
% folders to process
% varargin: 'standard', 'global_remap', 'rate_remap'
% outputs plot as well

parameters= list_of_parameters;
if isempty(varargin)
    varargin{1}= 'standard';
end

data_folder= pwd;
disp('processing folders for ...')
disp(varargin)
decoding_error= table;
counts= 1;
for this_folder=1:length(folders)
    disp(folders{this_folder})
    for this_method= 1:length(varargin)
        switch varargin{this_method}
            case 'standard'  
                cd([data_folder '\' folders{this_folder}])
                type= {'standard'};
            case 'global_remap'
                cd(['CONTROLS\global_remapped\' folders{this_folder}]); 
                type= {'global_remap'};
            case 'rate_remap'
                cd(['CONTROLS\rate_remapped\' folders{this_folder}]);   
                type= {'rate_remap'};   
        end
        load('estimated_position_leave_one_out.mat');
        if ~isfield(estimated_position_test,'percent_small_error') 
            if~contains(type{1},'standard')
            bayesian_decoding_error('method','leave one out','control',type{1});
            else
            bayesian_decoding_error('method','leave one out');
            end
            load('estimated_position_leave_one_out.mat');
        end
        if exist('estimated_position.mat')~=2
            spike_count([],[],[],'Y');
            bayesian_decoding([],[],'Y');
        end
        load('estimated_position.mat');
        num_tracks= length(estimated_position_test);
        decoding_error.session(counts:counts+num_tracks-1)= repmat(folders(this_folder),num_tracks,1);
        decoding_error.track(counts:counts+num_tracks-1)= [1:num_tracks];
        decoding_error.type(counts:counts+num_tracks-1)= type;
        % use function to get variables of interest
        [median_error,percent_small,accuracy]= get_decoding_measures(estimated_position_test,estimated_position,parameters);
        decoding_error.median_error(counts:counts+num_tracks-1)= median_error;    
        decoding_error.percent_small(counts:counts+num_tracks-1)= percent_small;    
        decoding_error.accuracy(counts:counts+num_tracks-1)= accuracy;          
        
        cd(data_folder);
        counts= counts+num_tracks;
    end
    
end

save('Tables\decoding_error_comparison.mat','decoding_error');

% figure('Color',[1 1 1]);
% groups=  [ones(size(decoding_error.standard)); 2*ones(size(decoding_error.rate_remapped)) ; 3*ones(size(decoding_error.global_remapped))];
% boxplot([decoding_error.standard; decoding_error.rate_remapped; decoding_error.global_remapped],groups,'PlotStyle','compact',...
%     'labels',{'standard','rate remapped','global remapped'},'LabelOrientation','horizontal','MedianStyle','line');
% a = get(get(gca,'children'),'children');
% set(a(4),'Color',[0.5 0.5 0.5]);
% set(a(5),'Color',[0.5 0.5 0.5]);
% set(a(6),'Color',[0.5 0.5 0.5]);
% set(a(7),'Color',parameters.colors.orange);
% set(a(8),'Color',parameters.colors.purple);
% set(a(9),'Color',parameters.colors.brick);
% set(a(7),'linewidth',20);
% set(a(8),'linewidth',20);
% set(a(9),'linewidth',20);
% hold on;
% x_pos= [1.2*ones(size(decoding_error.standard)); 2.2*ones(size(decoding_error.rate_remapped)) ; 3.2*ones(size(decoding_error.global_remapped))];
% x_pos= x_pos + 0.2*rand(length(x_pos),1);
% plot(x_pos, [decoding_error.standard; decoding_error.rate_remapped; decoding_error.global_remapped],'k.','MarkerSize',10);
% ylabel('median decoding error');
% [~,~,stats] = anova1([decoding_error.standard; decoding_error.rate_remapped; decoding_error.global_remapped],groups,'off');
% [c,~,~,~] = multcompare(stats,'Display','off'); % if you want to check exact values
% y_max= max([decoding_error.standard; decoding_error.rate_remapped]);
% plot([1 2],[y_max+5 y_max+5],'k','LineWidth',1.5);
% text(1.4,y_max+8,'n.s')
% y_max2= max(decoding_error.global_remapped);
% plot([1.5 3],[y_max2+5 y_max2+5],'k','LineWidth',1.5);
% text(2.2,y_max2+7,'***')
% ylim([0 y_max2+15])
end

function [median_error,percent_small,accuracy]= get_decoding_measures(est_wt,estimated_position,parameters)
    
        median_error= [est_wt.median_error]';
        percent_small= [est_wt.percent_small_error]';

        % re get number classification accuracy for GLMs as well
        load('extracted_position.mat');
        all_prob= vertcat(estimated_position.run); % look at all probabilities
        bayes_bias= vertcat(estimated_position.run_bias);
        [max_prob,max_prob_pos]= max(all_prob,[],1); % get index, first track will be first rows etc.. 
        no_spk_bins= find(sum(all_prob,1)==0);
        max_prob_pos(no_spk_bins)=[];
        max_prob(no_spk_bins)= [];
        
        num_bins_tracks= cellfun(@(x) size(x,1), {estimated_position.run});
        delim_bins_tracks(:,1)= [1 cumsum(num_bins_tracks(1:end-1))+1];
        delim_bins_tracks(:,2)= cumsum(num_bins_tracks);

        track_classification = NaN(size(max_prob_pos));
        for track_id=1:length(estimated_position)
            % classify which track max prob falls in
            indices= max_prob_pos >= delim_bins_tracks(track_id,1) & max_prob_pos <= delim_bins_tracks(track_id,2);
            track_classification(indices)= track_id;
        end

        for track_id=1:length(estimated_position)
            track_classification_tmp= track_classification;
            t_bins= estimated_position(track_id).run_time_centered;
            t_bins(no_spk_bins)=[];
            % to be fair, only take times where v_cm > 5cm/s
%             estimated_position(track_id).run_time_centered= [bayesian_spike_count.run_time_centered];
            speed_discrete= interp1(position.t, position.v_cm, t_bins, 'nearest');
            pos_discrete= interp1(position.t, estimated_position(track_id).discrete_position, t_bins, 'nearest');
            low_speed_idx= speed_discrete < parameters.speed_threshold_laps;
            track_classification_tmp(low_speed_idx)= NaN;
            pos_discrete(low_speed_idx)= NaN;

            accuracy(track_id)= numel(find(track_classification_tmp(~isnan(pos_discrete)) == track_id))./numel(track_classification_tmp(~isnan(pos_discrete)));
        end

end