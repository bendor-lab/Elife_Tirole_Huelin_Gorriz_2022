

function plot_control_track_peakFR(epochs,threshold_direction)
% Threshold direction (above/below [str]): indicates whether the selected cells are below or
    % above the threshold. 

if isempty(epochs)
    epochs = {'PRE','POST'};
end

load('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\rate_remapping_analysis_TRACK_PAIRS_wcorr.mat')
load('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\Tables\subsets_of_cells.mat');

thresholds = 1:1:10;

if length(find(contains(subset_of_cells.subset,['Track peakFR ' threshold_direction]))) < length(thresholds) %if subsets don't exist

    % Prepare cells across sessions
    all_cells = [remapping_raw(:).new_ID];
    [unique_cells,ia,~] = unique(all_cells,'stable');
    all_T1_peakFr = [remapping_raw(:).raw_peak_BAYESIAN_plfield_1];
    all_T1_peakFr = all_T1_peakFr(ia);
    all_T2_peakFr = [remapping_raw(:).raw_peak_BAYESIAN_plfield_2];
    all_T2_peakFr = all_T2_peakFr(ia);    
    
    for t = 1 : length(thresholds)

        if strcmp(threshold_direction,'below')
            % For each track, find cells that have peak FR below threshold
            t1_cell_IDs = all_cells(all_T1_peakFr < thresholds(t));
            t2_cell_IDs = all_cells(all_T2_peakFr < thresholds(t));
        elseif strcmp(threshold_direction,'above')
            % For each track, find cells that have peak FR above threshold
            t1_cell_IDs = all_cells(all_T1_peakFr > thresholds(t));
            t2_cell_IDs = all_cells(all_T2_peakFr > thresholds(t));
        end

        subset_cell_IDs = unique([t1_cell_IDs t2_cell_IDs]); % get cell IDs
        name_subset = ['Track peakFR ' threshold_direction ' ' num2str(thresholds(t)) ' spks/s'];
        % save in subset_of_cells a vector with cell IDs of cells passing the
        % threshold
        if any(strcmp(subset_of_cells.subset,name_subset)) %if already exists
            row_id = find(strcmp(subset_of_cells.subset,name_subset));
            subset_of_cells.cell_IDs{row_id} = subset_cell_IDs;
        else
            new_row = size(subset_of_cells,1)+1;
            subset_of_cells.subset{new_row} = name_subset;
            subset_of_cells.cell_IDs{new_row} = subset_cell_IDs;
        end

    end
    save('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\Tables\subsets_of_cells.mat','subset_of_cells')

end



% Get F and P-values for each subset
for t = 1 : length(thresholds)
    name_subset = ['Track peakFR ' threshold_direction ' ' num2str(thresholds(t)) ' spks/s'];
    if isempty(subset_of_cells.cell_IDs{ismember(subset_of_cells.subset,name_subset)})
        [pval(:,t),Fstat(:,t)] = deal(NaN);
        s_intcp(:,t,:) = [nan;nan];
        number_of_cells(t) = 0;
        continue
    end
    [pval(:,t),Fstat(:,t), s_intcp(:,t,:)] = plot_rate_remapping_NEW('x_var',{'place_field_diff'},'y_var',{'mean_max_FR_replay_diff'},...
        'subset',{'stable cells laps',name_subset},'y_label','Peak replay rate diff (Hz)','x_label','Peak in-field rate diff(Hz)','epochs',{epochs});
    close all
    number_of_cells(t) = length(intersect(subset_of_cells.cell_IDs{strcmp(subset_of_cells.subset,'stable cells laps')},...
        subset_of_cells.cell_IDs{strcmp(subset_of_cells.subset,name_subset)}));
end
%slopes = s_intcp(:,:,1);
%intercepts = s_intcp(:,:,2);

% Plot distribution
% f1 = figure('units','normalized','Color','w');
% for k = 1 : size(pval,1)
%     ax(k)= subplot(1,length(epochs),k);
%     plot(Fstat(k,:),'k.-','LineWidth',1.5,'MarkerSize',20);
%     ylabel('F Statistic')
%     hold on;
%     xlabel('Peak FR threshold (spks/s)')
%     yyaxis right
%     plot(log10(pval(k,:)),'r-','LineWidth',1.5,'MarkerSize',20);
%     plot(slopes(k,:).*sum(gausswin(0.25/(1/1000))),'Color',[0.5 0.5 0.5],'LineWidth',1.5,'MarkerSize',20,'Marker','.');
%     legend({'FStat','p value','slope'},'Box','off','Location','southeast')
%     ylabel('log_1_0(pval), slope')
% end
% 
% run_format_settings(f1,'match_ax')
% 
% f2 = figure;
% plot(number_of_cells,'k')
% run_format_settings(f2,'match_ax')
% 
% % Plot distribution
% f3 = figure('units','normalized','Color','w');
% for k = 1 : size(pval,1)
%     ax(k)= subplot(length(epochs),1,k);
%     p1 = plot(Fstat(k,:),'k.-','LineWidth',1.5,'MarkerSize',20);
%     ylabel('F Statistic')
%     hold on;
%     xlabel('Peak FR threshold (spks/s)')
%     yyaxis right
%     p2 = plot(log10(pval(k,:)),'Color',[0 .5 1],'LineWidth',1.5,'MarkerSize',20);
%     ylabel('log_1_0(pval)')
%     plot([min(xlim(ax(k))) max(xlim(ax(k)))],[log(0.05) log(0.05)],':','Color',[0.6 0.6 0.6],'LineWidth',1.5)    
%     legend([p1,p2],{'FStat','p value'},'Box','off','Location','southeast')
% 
% end
% 
% run_format_settings(f3,'match_ax')

% Plot distribution

% logged_pval = log10(pval);
% step_size = 1;
% min_val = min(min(logged_pval));
% max_val = max(max(logged_pval));
% new_lab = [min_val log10(0.05) max_val];% change colormap vals to log
% if max_val < log10(.05)
%     max_val = log10(.05);
%    % new_lab = [min_val ((max_val-min_val)/2)+min_val max_val];% change colormap vals to log
%         new_lab = [min_val log10(0.025) max_val];% change colormap vals to log
% 
% end  
% steps = min_val:step_size:max_val;
% steps(end+1) = steps(end) + step_size;
% temp = reshape(logged_pval',[(size(logged_pval,1)*size(logged_pval,2)),1]);
% [~,~,bins] = histcounts(temp,steps);
% colmap = zeros(length(steps),3);
% colmap(:,1) = linspace(1,0,length(steps)); %0,1
% colmap(:,2) = linspace(0,0,length(steps));
% colmap(:,3) = linspace(0,1,length(steps)); %1,0
% cdata_diff = zeros(length(bins),3);
% for k = 1 : length(steps)
%     cdata_diff(bins == k,:) = repmat(colmap(k,:),length(find(bins == k)==1),1);
% end
% 

% create colormap
steps = 0:0.001:0.05;
colmap = zeros(length(steps),3);
colmap(:,1) = linspace(1,0.7,length(steps)); 
colmap(:,2) = linspace(0,0,length(steps));
colmap(:,3) = linspace(0,0.7,length(steps));
added_steps= linspace(0.06,1,5);
colmap2 = zeros(length(added_steps),3);
colmap2(:,1)= linspace(0.5,0,length(added_steps));
colmap2(:,2) = linspace(0,0,length(added_steps));
colmap2(:,3) = linspace(0.5,1,length(added_steps));

colmap= ([colmap ; colmap2]);
new_steps= [steps added_steps];
[~,~,bin_id]= histcounts(pval',new_steps);
if any(ismember(bin_id,0))
    pval_colmap = colmap(bin_id(~ismember(bin_id,0)),:);
    nan_idx = find(ismember(bin_id,0));
    if nan_idx == 1
        pval_colmap = [0.3 0.3 0.3; pval_colmap];
    elseif nan_idx == length(pval)
        pval_colmap = [pval_colmap; 0.3 0.3 0.3];
    else
        pre_nan_idcs = 1: nan_idx-1;
        post_nan_idcs =  nan_idx+1:length(pval);
        pval_colmap = [pval_colmap(pre_nan_idcs,:); 0.3 0.3 0.3;pval_colmap(post_nan_idcs,:)];
    end    
else
    pval_colmap = colmap(bin_id,:);
end

f3 = figure('units','normalized','Color','w');
logged_pval = pval;
cdata_diff = pval_colmap;
for k = 1 : size(logged_pval,1)
    if k == 1
        this_col = cdata_diff(1:size(logged_pval,2),:);
    else
        new_s = size(logged_pval,2) * (k-1);
        this_col = cdata_diff(new_s+1:new_s+size(logged_pval,2),:);
    end
    %ax(k)= subplot(1,length(epochs),k);
    p1 = plot(Fstat(k,:),'k','LineWidth',1.5);
    hold on
    scatter(1:length(Fstat(k,:)),Fstat(k,:),50,this_col,'filled');
    ylabel('F Statistic')
    ytickformat('%i')
    yyaxis right
    p2 = plot(1:length(number_of_cells),number_of_cells,'Color', [0.4 0.4 0.4],'LineWidth',1.5);
    set(gca,{'ycolor'},{[0.4 .4 .4]})
    ylabel('number of cells','Rotation',270)
    if strcmp(threshold_direction,'above')
        xlabel({'Min peak in-field firing rate';'threshold (spks/s)'})
    else
        xlabel({'Max peak in-field firing rate';'threshold (spks/s)'})
    end
   % plot([min(xlim(ax(k))) max(xlim(ax(k)))],[log(0.05) log(0.05)],':','Color',[0.6 0.6 0.6],'LineWidth',1.5)    
    legend([p1,p2],{'FStat','number of cells'},'Box','off','Location','northeast')

end
ax=gca;
if strcmp(threshold_direction,'above')
    set(ax,'XLim',[0 10.5])
    yyaxis left
    set(ax,'YLim',[20 max(ylim)])
else
    set(ax,'XLim',[1 10.5])
    set(ax,'YLim',[0 max(ylim)+20])
end
h = colorbar;
colormap(colmap)
run_format_settings(f3,'match_ax')
%h.TickLabels = num2cell(round(new_lab,1));
h.Ticks= [0 find(new_steps==0.05)/size(colmap,1) 1]; 
h.TickLabels = [0 0.05 1];
ylabel(h,'pval','Rotation',270);
end