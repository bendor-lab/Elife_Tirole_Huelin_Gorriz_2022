function plot_cumulative_replay
load replay_CUMULATIVE
parameters=list_of_parameters;
number_of_tracks=length(out);
figure
for track=1:number_of_tracks
subplot(4,1,1)
hold on
plot((out(track).sleepPRE_time_bins(1:end-1)+out(track).bin_width/2)/60,out(track).sleepPRE_HIST,parameters.plot_color_line{track});
xlabel('cumulative sleep time (min)')
title('Replay during PRE')
ylabel('number of events / min') 
subplot(4,1,2)
hold on
plot((out(track).sleepPOST_time_bins(1:end-1)+out(track).bin_width/2)/60,out(track).sleepPOST_HIST,parameters.plot_color_line{track});
xlabel('cumulative sleep time (min)')
title('Replay during POST')
ylabel('number of events / min') 
subplot(4,1,3)
hold on
plot((out(track).sleepPOST_time_bins(1:end-1)+out(track).bin_width/2)/60,out(track).sleepPOST_CUMULATIVE,parameters.plot_color_line{track});
xlabel('cumulative sleep time (min)')
title('cumulative Replay during POST')
ylabel('total number of events/ session')
end
subplot(4,1,4)
if number_of_tracks==4
    plot((out(track).sleepPOST_time_bins(1:end-1)+out(track).bin_width/2)/60,out(1).sleepPOST_CUMULATIVE+out(3).sleepPOST_CUMULATIVE-out(2).sleepPOST_CUMULATIVE-out(4).sleepPOST_CUMULATIVE,'b');
ylabel('replay bias first track(+),second track(-)')
else
    for track=1:number_of_tracks
        sleepPOST_CUMULATIVE(track,:)=out(track).sleepPOST_CUMULATIVE;
    end
    for track=1:number_of_tracks
        plot((out(track).sleepPOST_time_bins(1:end-1)+out(track).bin_width/2)/60,sleepPOST_CUMULATIVE(track,:)-mean(sleepPOST_CUMULATIVE,1),parameters.plot_color_line{track});
    end
    ylabel('replay bias , current track vs. mean')
end
xlabel('cumulative sleep time (min)')

%plot a line for the duration of sleep in each subject...
subplot(4,1,1)
hold on
for i=1:length(out(1).sleepPRE_time_limit)
plot([0 out(1).sleepPRE_time_limit(i)]/60,[i i]/10,'k')
end
subplot(4,1,2)
hold on
for i=1:length(out(1).sleepPOST_time_limit)
plot([0 out(1).sleepPOST_time_limit(i)]/60,[i i]/10,'k')
end
subplot(4,1,3)
hold on
for i=1:length(out(1).sleepPOST_time_limit)
plot([0 out(1).sleepPOST_time_limit(i)]/60,[i i]/10,'k')
end
subplot(4,1,4)
hold on
plot([0 out(1).sleepPOST_time_limit(i)]/60,[0 0],'k--')
end
