function plot_all_tracks(varargin)

cd 'X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data'
load folders_to_process_remapping

parameters = list_of_parameters;

plfields1_PRE = [];
plfields2_PRE = [];
plfields1_POST =[];
plfields2_POST = [];
mean_spikes_PRE_1 = [];
mean_spikes_nonzero_PRE_1 = [];
mean_spikes_PRE_2 = [];
mean_spikes_nonzero_PRE_2 = [];
mean_spikes_POST_1 = [];
mean_spikes_nonzero_POST_1 = [];
mean_spikes_nonzero_POST_2 =[]
mean_spikes_POST_2 = [];
f1= figure;
f2=figure;
for i = 1 : length(folders)
    
    cd(['X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\CONTROLS\Tracks individually\' folders{i}])
    
    % Go back to original data folder to find cells used in common cells analysis
    curr_folder =  pwd;
    foldername = strsplit(pwd,'\');
    cd(['X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\' foldername{end}])
    if ~isempty(varargin)
        switch varargin{1}
            case 'wcorr'
                load('rate_remapping_analysis_TRACK_PAIRS_wcorr');
            case 'spearman'
                load('rate_remapping_analysis_TRACK_PAIRS_spearman')
        end
    end
    t2 = remapping;
    clear remapping
    
    cd(curr_folder)
    
    % Load track data
    if ~isempty(varargin)
        switch varargin{1}
            case 'wcorr'
                load  all_cells_analysis_single_track_wcorr
            case 'spearman'
                load all_cells_analysis_single_track_spearman
        end
    end
    
    if ~isempty(remapping(1).epoch(1).experiment)
        plfields1_PRE = [plfields1_PRE remapping(1).epoch(1).raw_peak_BAYESIAN_plfield_1];
        plfields2_PRE = [plfields2_PRE remapping(2).epoch(1).raw_peak_BAYESIAN_plfield_1];
        mean_spikes_PRE_1 = [mean_spikes_PRE_1; remapping(1).epoch(1).track1_mean_replay_spikes];
        mean_spikes_nonzero_PRE_1 = [mean_spikes_nonzero_PRE_1; remapping(1).epoch(1).track1_mean_replay_spikes_nonZero];
        mean_spikes_PRE_2 = [mean_spikes_PRE_2; remapping(2).epoch(1).track1_mean_replay_spikes];
        mean_spikes_nonzero_PRE_2 = [mean_spikes_nonzero_PRE_2; remapping(2).epoch(1).track1_mean_replay_spikes_nonZero];
    end
    
    plfields1_POST = [plfields1_POST remapping(1).epoch(2).raw_peak_BAYESIAN_plfield_1];
    plfields2_POST = [plfields2_POST remapping(2).epoch(2).raw_peak_BAYESIAN_plfield_1];
    mean_spikes_POST_1 = [mean_spikes_POST_1; remapping(1).epoch(2).track1_mean_replay_spikes];
    mean_spikes_nonzero_POST_1 = [mean_spikes_nonzero_POST_1; remapping(1).epoch(2).track1_mean_replay_spikes_nonZero];
    mean_spikes_POST_2 = [mean_spikes_POST_2; remapping(2).epoch(2).track1_mean_replay_spikes];
    mean_spikes_nonzero_POST_2 = [mean_spikes_nonzero_POST_2; remapping(2).epoch(2).track1_mean_replay_spikes_nonZero];
    
    
    
    %%%%% PLOT T1
    figure(f1)
    
    % PRE
    subplot(2,2,1) 
    % Find cells used in common cells analysis
    common_cells = t2(1).ID_active_cells_during_replay;
    [~,diff_analysis_idx,~] = intersect(remapping(1).epoch(1).ID_active_cells_during_replay,common_cells);
    small_FR_change = find(t2(1).place_field_diff < 1);
    [~,FR_change_idx,~] = intersect(remapping(1).epoch(1).ID_active_cells_during_replay,t2(1).ID_active_cells_during_replay(small_FR_change));
    
    hold on
    plot(remapping(1).epoch(1).raw_peak_BAYESIAN_plfield_1, remapping(1).epoch(1).track1_mean_replay_spikes_nonZero,'o','MarkerEdgeColor',[0.2 0.4 0.8]);
    hold on
    plot(remapping(1).epoch(1).raw_peak_BAYESIAN_plfield_1(diff_analysis_idx), remapping(1).epoch(1).track1_mean_replay_spikes_nonZero(diff_analysis_idx),'o','MarkerEdgeColor',[0.2 0.4 0.8],'MarkerFaceColor',[0.2 0.4 0.8]);
    plot(remapping(1).epoch(1).raw_peak_BAYESIAN_plfield_1(FR_change_idx), remapping(1).epoch(1).track1_mean_replay_spikes_nonZero(FR_change_idx),'o','MarkerEdgeColor','m','MarkerFaceColor',[0.2 0.4 0.8]);
    if i == 5
        lm = fitlm(plfields1_PRE,mean_spikes_nonzero_PRE_1,'linear')
        [p,~,~] = coefTest(lm)
        linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(1,1) = p;
        x=[min(plfields1_PRE) max(plfields1_PRE)];
        b=lm.Coefficients.Estimate';
        hold on
        plot(x,polyval(fliplr(b),x),'k--');
        title(num2str(p))
    end
    
    subplot(2,2,3)
    hold on
    plot(remapping(1).epoch(1).raw_peak_BAYESIAN_plfield_1, remapping(1).epoch(1).track1_mean_replay_spikes,'o','MarkerEdgeColor',[0.2 0.4 0.8]);
    hold on
    plot(remapping(1).epoch(1).raw_peak_BAYESIAN_plfield_1(diff_analysis_idx), remapping(1).epoch(1).track1_mean_replay_spikes(diff_analysis_idx),'o','MarkerEdgeColor',[0.2 0.4 0.8],'MarkerFaceColor',[0.2 0.4 0.8]);
    plot(remapping(1).epoch(1).raw_peak_BAYESIAN_plfield_1(FR_change_idx), remapping(1).epoch(1).track1_mean_replay_spikes(FR_change_idx),'o','MarkerEdgeColor','m','MarkerFaceColor',[0.2 0.4 0.8]);
    if i == 5
        lm = fitlm(plfields1_PRE,mean_spikes_PRE_1,'linear')
        [p,~,~] = coefTest(lm)
        linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(1,3) = p;
        x=[min(plfields1_PRE) max(plfields1_PRE)];
        b=lm.Coefficients.Estimate';
        hold on
        plot(x,polyval(fliplr(b),x),'k--');
        title(num2str(p))
    end
    
    % POST
    subplot(2,2,2)
    % Find cells used in common cells analysis
    common_cells = t2(2).ID_active_cells_during_replay;
    [~,diff_analysis_idx,~] = intersect(remapping(1).epoch(2).ID_active_cells_during_replay,common_cells);
    small_FR_change = find(t2(2).place_field_diff < 1);
    [~,FR_change_idx,~] = intersect(remapping(1).epoch(2).ID_active_cells_during_replay,t2(2).ID_active_cells_during_replay(small_FR_change));

    hold on
    plot(remapping(1).epoch(2).raw_peak_BAYESIAN_plfield_1,remapping(1).epoch(2).track1_mean_replay_spikes_nonZero,'o','MarkerEdgeColor',[0.2 0.4 0.8]);
    hold on
    plot(remapping(1).epoch(2).raw_peak_BAYESIAN_plfield_1(diff_analysis_idx), remapping(1).epoch(2).track1_mean_replay_spikes_nonZero(diff_analysis_idx),'o','MarkerEdgeColor',[0.2 0.4 0.8],'MarkerFaceColor',[0.2 0.4 0.8]);
    plot(remapping(1).epoch(2).raw_peak_BAYESIAN_plfield_1(small_FR_change), remapping(1).epoch(2).track1_mean_replay_spikes_nonZero(small_FR_change),'o','MarkerEdgeColor','m','MarkerFaceColor',[0.2 0.4 0.8]);
    if i == 5
        lm = fitlm(plfields1_POST,mean_spikes_nonzero_POST_1,'linear')
        [p,~,~] = coefTest(lm)
        linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(1,2) = p;
        x=[min(plfields1_POST) max(plfields1_POST)];
        b=lm.Coefficients.Estimate';
        hold on
        plot(x,polyval(fliplr(b),x),'k--');
        title(num2str(p))
    end
    
    subplot(2,2,4)
    hold on
    plot(remapping(1).epoch(2).raw_peak_BAYESIAN_plfield_1, remapping(1).epoch(2).track1_mean_replay_spikes,'o','MarkerEdgeColor',[0.2 0.4 0.8]);
    hold on
    plot(remapping(1).epoch(2).raw_peak_BAYESIAN_plfield_1(diff_analysis_idx), remapping(1).epoch(2).track1_mean_replay_spikes(diff_analysis_idx),'o','MarkerEdgeColor',[0.2 0.4 0.8],'MarkerFaceColor',[0.2 0.4 0.8]);
    plot(remapping(1).epoch(2).raw_peak_BAYESIAN_plfield_1(small_FR_change), remapping(1).epoch(2).track1_mean_replay_spikes(small_FR_change),'o','MarkerEdgeColor','m','MarkerFaceColor',[0.2 0.4 0.8]);
    if i == 5
        lm = fitlm(plfields1_POST,mean_spikes_POST_1,'linear')
        [p,~,~] = coefTest(lm)
        linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(1,4) = p;
        x=[min(plfields1_POST) max(plfields1_POST)];
        b=lm.Coefficients.Estimate';
        hold on
        plot(x,polyval(fliplr(b),x),'k--');
        title(num2str(p))
    end
    
    %%%% PLOT T2
    
    figure(f2)
    subplot(2,2,1)
    % Find cells used in common cells analysis
    common_cells = t2(1).ID_active_cells_during_replay;
    [~,diff_analysis_idx,~] = intersect(remapping(2).epoch(1).ID_active_cells_during_replay,common_cells);
    small_FR_change = find(t2(1).place_field_diff < 1);
    [~,FR_change_idx,~] = intersect(remapping(2).epoch(1).ID_active_cells_during_replay,t2(1).ID_active_cells_during_replay(small_FR_change));

    hold on
    plot(remapping(2).epoch(1).raw_peak_BAYESIAN_plfield_1, remapping(2).epoch(1).track1_mean_replay_spikes_nonZero,'o','MarkerEdgeColor',[0.2 0.4 0.8]);
    hold on
    plot(remapping(2).epoch(1).raw_peak_BAYESIAN_plfield_1(diff_analysis_idx), remapping(2).epoch(1).track1_mean_replay_spikes_nonZero(diff_analysis_idx),'o','MarkerEdgeColor',[0.2 0.4 0.8],'MarkerFaceColor',[0.2 0.4 0.8]);
    plot(remapping(2).epoch(1).raw_peak_BAYESIAN_plfield_1(FR_change_idx), remapping(2).epoch(1).track1_mean_replay_spikes_nonZero(FR_change_idx),'o','MarkerEdgeColor','m','MarkerFaceColor',[0.2 0.4 0.8]);
    if i == 5
        lm = fitlm(plfields2_PRE,mean_spikes_nonzero_PRE_2,'linear')
        [p,~,~] = coefTest(lm);
        linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(2,1) = p;
        x=[min(plfields2_PRE) max(plfields2_PRE)];
        b=lm.Coefficients.Estimate';
        hold on
        plot(x,polyval(fliplr(b),x),'k--');
        title(num2str(p))
    end
    
    subplot(2,2,3)
    hold on
    plot(remapping(2).epoch(1).raw_peak_BAYESIAN_plfield_1, remapping(2).epoch(1).track1_mean_replay_spikes,'o','MarkerEdgeColor',[0.2 0.4 0.8]);
    hold on
    plot(remapping(2).epoch(1).raw_peak_BAYESIAN_plfield_1(diff_analysis_idx), remapping(2).epoch(1).track1_mean_replay_spikes(diff_analysis_idx),'o','MarkerEdgeColor',[0.2 0.4 0.8],'MarkerFaceColor',[0.2 0.4 0.8]);
    plot(remapping(2).epoch(1).raw_peak_BAYESIAN_plfield_1(FR_change_idx), remapping(2).epoch(1).track1_mean_replay_spikes(FR_change_idx),'o','MarkerEdgeColor','m','MarkerFaceColor',[0.2 0.4 0.8]);
    if i == 5
        lm = fitlm(plfields2_PRE,mean_spikes_PRE_2,'linear')
        [p,~,~] = coefTest(lm);
        linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(2,3) = p;
        x=[min(plfields2_PRE) max(plfields2_PRE)];
        b=lm.Coefficients.Estimate';
        hold on
        plot(x,polyval(fliplr(b),x),'k--');
        title(num2str(p))
    end
    
    %POST
    subplot(2,2,2)
    % Find cells used in common cells analysis
    common_cells = t2(2).ID_active_cells_during_replay;
    [~,diff_analysis_idx,~] = intersect(remapping(2).epoch(2).ID_active_cells_during_replay,common_cells);
     small_FR_change = find(t2(2).place_field_diff < 1);
    [~,FR_change_idx,~] = intersect(remapping(2).epoch(2).ID_active_cells_during_replay,t2(2).ID_active_cells_during_replay(small_FR_change));

    hold on
    plot(remapping(2).epoch(2).raw_peak_BAYESIAN_plfield_1,remapping(2).epoch(2).track1_mean_replay_spikes_nonZero,'o','MarkerEdgeColor',[0.2 0.4 0.8]);
    hold on
    plot(remapping(2).epoch(2).raw_peak_BAYESIAN_plfield_1(diff_analysis_idx), remapping(2).epoch(2).track1_mean_replay_spikes_nonZero(diff_analysis_idx),'o','MarkerEdgeColor',[0.2 0.4 0.8],'MarkerFaceColor',[0.2 0.4 0.8]);
    plot(remapping(2).epoch(2).raw_peak_BAYESIAN_plfield_1(FR_change_idx), remapping(2).epoch(2).track1_mean_replay_spikes_nonZero(FR_change_idx),'o','MarkerEdgeColor','m','MarkerFaceColor',[0.2 0.4 0.8]);
    if i == 5
        lm = fitlm(plfields2_POST,mean_spikes_nonzero_POST_2,'linear')
        linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(2,2) = p;
        [p,~,~] = coefTest(lm);
        x=[min(plfields2_POST) max(plfields2_POST)];
        b=lm.Coefficients.Estimate';
        hold on
        plot(x,polyval(fliplr(b),x),'k--');
        title(num2str(p))
    end
    
    subplot(2,2,4)
    hold on
    plot(remapping(2).epoch(2).raw_peak_BAYESIAN_plfield_1, remapping(2).epoch(2).track1_mean_replay_spikes,'o','MarkerEdgeColor',[0.2 0.4 0.8]);
    hold on
    plot(remapping(2).epoch(2).raw_peak_BAYESIAN_plfield_1(diff_analysis_idx), remapping(2).epoch(2).track1_mean_replay_spikes(diff_analysis_idx),'o','MarkerEdgeColor',[0.2 0.4 0.8],'MarkerFaceColor',[0.2 0.4 0.8]);
    plot(remapping(2).epoch(2).raw_peak_BAYESIAN_plfield_1(FR_change_idx), remapping(2).epoch(2).track1_mean_replay_spikes(FR_change_idx),'o','MarkerEdgeColor','m','MarkerFaceColor',[0.2 0.4 0.8]);
    if i == 5
        lm = fitlm(plfields2_POST,mean_spikes_POST_2,'linear')
        [p,~,~] = coefTest(lm);
        linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(2,4) = p;
        x=[min(plfields2_POST) max(plfields2_POST)];
        b=lm.Coefficients.Estimate';
        hold on
        plot(x,polyval(fliplr(b),x),'k--');
        title(num2str(p))
    end
    
end


    figure(f1)
    f1.Name = ['T1_' varargin{1}];
    ax1 = subplot(2,2,1)
    title(['sleep- PRE nonzero , pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(1,1))])
    ylabel({'Mean # spikes T1 replay'});
    xlabel('T1 place field peak FR')
    box off
    ax1.FontSize = 16;

    ax2 = subplot(2,2,2)
    title(['sleep- POST nonzero , pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(1,2))])
    ylabel({'Mean # spikes T1 replay'});
    xlabel('T1 place field peak FR')
    box off
    ax2.FontSize = 16;

    ax3 = subplot(2,2,3)
    title(['sleep- PRE, pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(1,3))])
    ylabel({'Mean # spikes T1 replay'});
    xlabel('T1 place field peak FR')
    box off
    ax3.FontSize = 16;

    ax4 = subplot(2,2,4)
    title(['sleep- POST, pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(1,4))])
    ylabel({'Mean # spikes T1 replay'});
    xlabel('T1 place field peak FR')
    box off
    ax4.FontSize = 16;

    figure(f2)
    f2.Name = ['T2_' varargin{1}];
    ax1 = subplot(2,2,1)
    title(['sleep- PRE nonzero , pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(2,1))])
    ylabel({'Mean # spikes T2 replay'});
    xlabel('T2 place field peak FR')
    box off
    ax1.FontSize = 16;

    ax2 = subplot(2,2,2)
    title(['sleep- POST nonzero , pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(2,2))])
    ylabel({'Mean # spikes T2 replay'});
    xlabel('T2 place field peak FR')
    box off
    ax2.FontSize = 16;

    ax3 = subplot(2,2,3)
    title(['sleep- PRE, pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(2,3))])
    ylabel({'Mean # spikes T2 replay'});
    xlabel('T2 place field peak FR')
    box off
    ax3.FontSize = 16;

    ax4 = subplot(2,2,4)
    title(['sleep- POST, pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(2,4))])
    ylabel({'Mean # spikes T2 replay'});
    xlabel('T2 place field peak FR')
    box off
    ax4.FontSize = 16;













end