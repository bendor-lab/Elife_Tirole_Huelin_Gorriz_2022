function out=cumulative_replay(folders)
if isempty(folders)
    a=pwd;
    folders={a};
    cd ..;
end
out.folders=folders;
cd(folders{1});
load sorted_replay

number_of_tracks=length(sorted_replay);
for track=1:number_of_tracks    
    out(track).experiment=[];
    out(track).sleepPRE=[];
    out(track).sleepPOST=[];
    out(track).sleepPRE_expt=[];
    out(track).sleepPOST_expt=[];
    out(track).sleepPRE_time_limit=[];
    out(track).sleepPOST_time_limit=[];
end
cd ..


for k=1:length(folders)
    cd(folders{k});
    load sorted_replay  
    load time_range
    for track=1:number_of_tracks
        out(track).sleepPRE=[out(track).sleepPRE sorted_replay(track).cumulative_event_time.sleepPRE];
        out(track).sleepPOST=[out(track).sleepPOST sorted_replay(track).cumulative_event_time.sleepPOST];
        out(track).sleepPRE_expt=[out(track).sleepPRE_expt k*ones(size(sorted_replay(track).cumulative_event_time.sleepPRE))];
        out(track).sleepPOST_expt=[out(track).sleepPOST_expt k*ones(size(sorted_replay(track).cumulative_event_time.sleepPOST))];
        out(track).sleepPRE_time_limit=[out(track).sleepPRE_time_limit max(max(time_range.sleepPRE_CUMULATIVE))];
        out(track).sleepPOST_time_limit=[out(track).sleepPOST_time_limit max(max(time_range.sleepPOST_CUMULATIVE))];
    end
    cd ..
end

bin_width=60;
for track=1:number_of_tracks
    out(track).bin_width=bin_width;
    out(track).sleepPRE(find(isnan(out(track).sleepPRE)))=[];
    out(track).sleepPOST(find(isnan(out(track).sleepPOST)))=[];
    
    out(track).sleepPRE_time_bins=0:out(track).bin_width:bin_width*ceil(max(out(track).sleepPRE_time_limit)/out(track).bin_width);
    if isempty(out(track).sleepPRE_time_bins)
        out(track).sleepPRE_HIST=[];
    else
        out(track).sleepPRE_HIST=histcounts(out(track).sleepPRE,out(track).sleepPRE_time_bins);
    end
    for i=1:length(out(track).sleepPRE_time_bins)-1
        out(track).number_of_sessionsPRE(i)=length(find(out(track).sleepPRE_time_limit-out(track).sleepPRE_time_bins(i)>=0)); %number of experiments where current time bin start is less than overall limit
    end
    out(track).sleepPOST_time_bins=0:bin_width:out(track).bin_width*ceil(max(out(track).sleepPOST_time_limit)/out(track).bin_width);
    if isempty(out(track).sleepPOST_time_bins)
        out(track).sleepPOST_HIST=[];
    else
        out(track).sleepPOST_HIST=histcounts(out(track).sleepPOST,out(track).sleepPOST_time_bins);
    end
    for i=1:length(out(track).sleepPOST_time_bins)-1
        out(track).number_of_sessionsPOST(i)=length(find(out(track).sleepPOST_time_limit-out(track).sleepPOST_time_bins(i)>=0)); %number of experiments where current time bin start is less than overall limit
    end
    if isempty(out(track).sleepPRE_HIST)
        out(track).sleepPRE_CUMULATIVE=[];
    else
        out(track).sleepPRE_CUMULATIVE=cumsum(out(track).sleepPRE_HIST./out(track).number_of_sessionsPRE);
        out(track).sleepPRE_HIST=smooth(out(track).sleepPRE_HIST./out(track).number_of_sessionsPRE,5);
    end
    if isempty(out(track).sleepPOST_HIST)
        out(track).sleepPOST_CUMULATIVE=[];
    else
        out(track).sleepPOST_CUMULATIVE=cumsum(out(track).sleepPOST_HIST./out(track).number_of_sessionsPOST);
        out(track).sleepPOST_HIST=smooth(out(track).sleepPOST_HIST./out(track).number_of_sessionsPOST,5);
    end
    
end


save replay_CUMULATIVE out

end