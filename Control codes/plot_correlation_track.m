
function plot_correlation_track(varargin)


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


parameters=list_of_parameters;

if ~isempty(varargin)
    switch varargin{1}
        case 'wcorr'
            load('all_cells_analysis_single_track_wcorr');
        case 'spearman'
            load('all_cells_analysis_single_track_spearman');
        otherwise
            load('all_cells_analysis_single_track');
    end
else
    load('rate_remapping_analysis_TRACK_PAIRS');
end

for t = 1 : size(remapping,2) % for each track
    foldername = strsplit(pwd,'\');
    f(t)= figure('units','normalized','outerposition',[0 0 1 1]);
    f(t).Name= strcat(varargin{1},'_',foldername{3},'_CONTROL_T',num2str(t));
    
    for epoch=1:size(remapping,2)
        data_across_tracks(epoch).place_field_peakFR=[];
        data_across_tracks(epoch).mean_spike=[];
        data_across_tracks(epoch).mean_spike_nonZero=[];
        data_across_tracks(epoch).track_id=[];
        [~,diff_analysis_idx,~] = intersect(remapping(t).epoch(epoch).ID_active_cells_during_replay,remapping(t).epoch(epoch).PRE_to_POST_active_cells);

        figure(f(t))
        subplot(size(remapping,2),2,epoch)
        plot(remapping(t).epoch(epoch).raw_peak_BAYESIAN_plfield_1, remapping(t).epoch(epoch).track1_mean_replay_spikes_nonZero,parameters.plot_color_symbol{1});
        hold on
        plot(remapping(t).epoch(epoch).raw_peak_BAYESIAN_plfield_1(diff_analysis_idx), remapping(t).epoch(epoch).track1_mean_replay_spikes_nonZero(diff_analysis_idx),'k*');
        subplot(size(remapping,2),2,epoch+2)
        plot(remapping(t).epoch(epoch).raw_peak_BAYESIAN_plfield_1, remapping(t).epoch(epoch).track1_mean_replay_spikes,parameters.plot_color_symbol{1});
        hold on
        plot(remapping(t).epoch(epoch).raw_peak_BAYESIAN_plfield_1(diff_analysis_idx), remapping(t).epoch(epoch).track1_mean_replay_spikes(diff_analysis_idx),'k*');
        
        index_non_NaNs=find(~isnan(remapping(t).epoch(epoch).raw_peak_BAYESIAN_plfield_1));
        data_across_tracks(epoch).place_field_peakFR=[data_across_tracks(epoch).place_field_peakFR; remapping(t).epoch(epoch).raw_peak_BAYESIAN_plfield_1(index_non_NaNs)];
        data_across_tracks(epoch).mean_spike=[data_across_tracks(epoch).mean_spike; remapping(t).epoch(epoch).track1_mean_replay_spikes(index_non_NaNs)];
        data_across_tracks(epoch).mean_spike_nonZero=[data_across_tracks(epoch).mean_spike_nonZero; remapping(t).epoch(epoch).track1_mean_replay_spikes_nonZero(index_non_NaNs)];
        
        
        if ~isempty(data_across_tracks(epoch).place_field_peakFR) & ~isempty(data_across_tracks(epoch).mean_spike)
            lm = fitlm(data_across_tracks(epoch).place_field_peakFR, data_across_tracks(epoch).mean_spike_nonZero,'linear')
            [p,~,~] = coefTest(lm);
            linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(epoch) = p;
            x=[min(data_across_tracks(epoch).place_field_peakFR) max(data_across_tracks(epoch).place_field_peakFR)];
            b=lm.Coefficients.Estimate';
            hold on
            subplot(size(remapping,2),2,epoch)
            plot(x,polyval(fliplr(b),x),'k--');
            
            
            
            lm = fitlm(data_across_tracks(epoch).place_field_peakFR, data_across_tracks(epoch).mean_spike,'linear')
            [p,~,~] = coefTest(lm);
            linear_fit_P_values_RATE_VS_PEAK_PF(epoch)= p;
            x=[min(data_across_tracks(epoch).place_field_peakFR) max(data_across_tracks(epoch).place_field_peakFR)];
            b=lm.Coefficients.Estimate';
            hold on
            subplot(size(remapping,2),2,epoch+2)
            plot(x,polyval(fliplr(b),x),'k--');
        end
        
    end
    
    
    subplot(2,2,1)
    title(['sleep- PRE nonzero , pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(1))])
    ylabel([{'Average spike number difference during replay between tracks'}; {'divided by number of events cell was active in)'}])
    xlabel('place field peak response difference between tracks')
    subplot(2,2,2)
    title(['sleep- POST nonzero , pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(2))])
    ylabel([{'Average spike number difference during replay between tracks'}; {'divided by number of events cell was active in)'}])
    xlabel('place field peak response difference between tracks')
    subplot(2,2,3)
    title(['sleep- PRE, pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF(1))])
    ylabel('Average spike number difference during replay between tracks')
    xlabel('place field peak response difference between tracks')
    subplot(2,2,4)
    title(['sleep- POST, pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF(2))])
    ylabel('Average spike number difference during replay between tracks')
    xlabel('place field peak response difference between tracks')
    
end

end