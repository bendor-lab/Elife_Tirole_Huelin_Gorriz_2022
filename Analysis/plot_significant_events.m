function plot_significant_events
load extracted_position
load extracted_sleep_state
load significant_replay_events


parameters=list_of_parameters;
number_of_tracks=length(significant_replay_events.track);
figure
ax(1)=subplot(6,1,1);
plot(position.t,position.v_cm);
xlabel('Time'); ylabel('velocity')

ax(2)=subplot(6,1,2);
sleep_periods = sleep_state.state_binned;
sleep_periods(sleep_periods==-1) = 0;
area(sleep_state.time,sleep_periods,'LineWidth',3,'FaceColor','g','FaceAlpha',0.2,'EdgeColor','g','EdgeAlpha',0.2);
hold on
plot(significant_replay_events.time_bin_centres,smooth(significant_replay_events.HIST,significant_replay_events.smoothing_bin_size),'k')
xlabel('Time'); ylabel('replay events')
legend('Sleep state','Amount of replay events')

ax(3)=subplot(6,1,3);
hold on
for track=1:number_of_tracks
    plot(significant_replay_events.time_bin_centres,smooth(significant_replay_events.track(track).HIST,significant_replay_events.smoothing_bin_size),parameters.plot_color_line{track})
end
plot(significant_replay_events.time_bin_centres,smooth(significant_replay_events.p_value_threshold*significant_replay_events.HIST,significant_replay_events.smoothing_bin_size),'k--')
plot(significant_replay_events.all_event_times(significant_replay_events.BAYESIAN_BIAS_excluded_events_index),...
    zeros(size(significant_replay_events.BAYESIAN_BIAS_excluded_events_index)),'ko')
xlabel('Time'); ylabel('number of replay events')
title('Weighted correlation - BAYESIAN BIAS METHOD')

%%%%%%OTHER OPTIONS%%%%%%%%%%%%%%
% plot(significant_replay_events.all_event_times(multi_tracks_index), zeros(size(significant_replay_events.all_event_times(multi_tracks_index))),'ko')
% plot(significant_replay_events.time_bin_centres,smooth(sig_event_info.track(track).HIST,significant_replay_events.smoothing_bin_size),parameters.plot_color_line{track})
% plot(significant_replay_events.time_bin_centres,smooth(sig_event_info.track(track).HIST_NO_MULTI,significant_replay_events.smoothing_bin_size),parameters.plot_color_line{track})
% plot(significant_replay_events.time_bin_centres,smooth(sig_event_info.track(track).HIST_MULTI,significant_replay_events.smoothing_bin_size),parameters.plot_color_dashed_line{track})
%
ax(4)=subplot(6,1,4);
hold on
for track=1:number_of_tracks
    plot(significant_replay_events.time_bin_centres_BIG(6:end-5),significant_replay_events.track(track).mean_best_score,parameters.plot_color_line{track})
end
xlabel('Time'); ylabel('mean score')
title('Replay score')

ax(5)=subplot(6,1,5);
hold on
for track=1:number_of_tracks
    plot(significant_replay_events.time_bin_centres_BIG(6:end-5),significant_replay_events.track(track).mean_bayesian_bias,parameters.plot_color_line{track})
end
xlabel('Time'); ylabel('mean bias')
title('bayesian bias')

ax(6)=subplot(6,1,6);
hold on
for track=1:number_of_tracks
    plot(significant_replay_events.time_bin_centres_BIG(6:end-5),significant_replay_events.track(track).min_log_pvalue,parameters.plot_color_line{track})
end
xlabel('Time'); ylabel('min log pvalue')
title('log p value')

linkaxes(ax,'x');

% figure
% for track=1:number_of_tracks
%     subplot(2,2,1)
%     hold on
%     plot(significant_replay_events.track(track).mean_best_score,significant_replay_events.track(track).median_log_pvalue,parameters.plot_color_symbol{track})
%     xlabel('mean best score'); ylabel('median log pvalue')
%     subplot(2,2,2)
%     hold on
%     plot(significant_replay_events.track(track).mean_best_score,significant_replay_events.track(track).mean_bayesian_bias,parameters.plot_color_symbol{track})
%     xlabel('mean best score'); ylabel('mean bayesian bias')
%     subplot(2,2,3)
%     hold on
%     plot(significant_replay_events.track(track).mean_bayesian_bias,significant_replay_events.track(track).median_log_pvalue,parameters.plot_color_symbol{track})
%     xlabel('mean bayesian bias'); ylabel('median log pvalue')
%     subplot(2,2,4)
%     hold on
% end