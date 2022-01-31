function [pvalue_versus_original pvalue_versus_shuffle] = DeLong_test(log_odd)
states = [-1,0,1,2];

track_id = log_odd.track;
epoch_track_index = [];
for k=1:length(states) % sort track id and behavioral states (pooled from 5 sessions)
    state_index = find(log_odd.behavioural_state==states(k));
    
    if k == 1 % PRE
        epoch_track_index{1}{1} = intersect(state_index,find(track_id==1));
        epoch_track_index{2}{1} = intersect(state_index,find(track_id==2));
    elseif k == 2 % POST
        epoch_track_index{1}{3} = intersect(state_index,find(track_id==1));
        epoch_track_index{2}{3} = intersect(state_index,find(track_id==2));
    elseif k == 3 % RUN Track 1
        epoch_track_index{1}{2} = intersect(state_index,find(track_id==1)); % Track 1 replay on Track 1 considered RUN Replay      
    elseif k == 4 % RUN Track 2
        epoch_track_index{2}{2} = intersect(state_index,find(track_id==2));        
    end
end

% Resampling so that no of T1 replay matches no of T2 replay 
for epoch=1:3 % sort track id and behavioral states
    difference = abs(length(epoch_track_index{1}{epoch}) - length(epoch_track_index{2}{epoch})); % find the difference between two tracks
    combined = [epoch_track_index{1}{epoch},epoch_track_index{2}{epoch}];
    
    if length(epoch_track_index{1}{epoch}) > length(epoch_track_index{2}{epoch})
        combined = [combined datasample(epoch_track_index{2}{epoch},difference)]; % combine and equate the sample size
        balanced_index{1}{epoch} = combined(1:length(combined)/2); % divide half into one
        balanced_index{2}{epoch} = combined(1+length(combined)/2:end); % divide another half into another one
    elseif length(epoch_track_index{1}{epoch}) < length(epoch_track_index{2}{epoch})
        if ~isempty(epoch_track_index{1}{epoch}) % If track 1 is empty then track 2 also empty
            combined = [datasample(epoch_track_index{1}{epoch},difference) combined];
            balanced_index{1}{epoch} = combined(1:length(combined)/2);
            balanced_index{2}{epoch} = combined(1+length(combined)/2:end);
        end
    end
end

for epoch = 1:3
    epoch_index{epoch} = [balanced_index{1}{epoch} balanced_index{2}{epoch}];
end

for epoch = 1:3
    % Original VS 20ms rate fixed
    sample = [];
    sample.ratings = [log_odd.normal_zscored.original(epoch_index{epoch}); log_odd.normal_zscored.rate_fixed(epoch_index{epoch})];
    sample.spsizes = [length(epoch_index{epoch})/2 length(epoch_index{epoch})/2];
    [aucs, delongcov] = fastDeLong(sample);
    z = diff(aucs) / sqrt(delongcov(1,1)+ delongcov(2,2)- 2*delongcov(1,2));
    pvalue_versus_original(epoch,1) =  2*normcdf(abs(z), 0, 1,'upper')
    
    % Original VS One bin
    sample = [];
    sample.ratings = [log_odd.normal_zscored.original(epoch_index{epoch}); log_odd.one_bin_zscored.original(epoch_index{epoch})];
    sample.spsizes = [length(epoch_index{epoch})/2 length(epoch_index{epoch})/2];
    [aucs, delongcov] = fastDeLong(sample);
    z = diff(aucs) / sqrt(delongcov(1,1)+ delongcov(2,2)- 2*delongcov(1,2));
    pvalue_versus_original(epoch,2) =  2*normcdf(abs(z), 0, 1,'upper')
    
    % Original Vs rate fixed global remapping
    sample = [];
    sample.ratings = [log_odd.normal_zscored.original(epoch_index{epoch}); log_odd.normal_zscored.rate_fixed_global_remapping(epoch_index{epoch})];
    sample.spsizes = [length(epoch_index{epoch})/2 length(epoch_index{epoch})/2];
    [aucs, delongcov] = fastDeLong(sample);
    z = diff(aucs) / sqrt(delongcov(1,1)+ delongcov(2,2)- 2*delongcov(1,2));
    pvalue_versus_original(epoch,3) =  2*normcdf(abs(z), 0, 1,'upper')
    % Original Vs one bin rate remapping
    sample = [];
    sample.ratings = [log_odd.normal_zscored.original(epoch_index{epoch}); log_odd.one_bin_zscored.rate_remapping(epoch_index{epoch})];
    sample.spsizes = [length(epoch_index{epoch})/2 length(epoch_index{epoch})/2];
    z = diff(aucs) / sqrt(delongcov(1,1)+ delongcov(2,2)- 2*delongcov(1,2));
    pvalue_versus_original(epoch,4) =  2*normcdf(abs(z), 0, 1,'upper')
    
%     sample = [];
%     sample.ratings = [log_odd.normal_zscored.original(epoch_index{epoch}); log_odd.one_bin_zscored.rate_fixed(epoch_index{epoch})];
%     sample.spsizes = [length(epoch_index{epoch})/2 length(epoch_index{epoch})/2];
%     z = abs(diff(aucs)) / sqrt(delongcov(1,1)+ delongcov(2,2)- 2*delongcov(1,2));
%     pvalue_versus_original(epoch,4) =  2*normcdf(abs(z), 0, 1,'upper')
    
    
    sample = [];
    
    % one bin rate remapping VS original
    sample.ratings = [log_odd.normal_zscored.original(epoch_index{epoch});log_odd.one_bin_zscored.rate_remapping(epoch_index{epoch})];
    sample.spsizes = [length(epoch_index{epoch})/2 length(epoch_index{epoch})/2];
    [aucs, delongcov] = fastDeLong(sample);
    z = diff(aucs) / sqrt(delongcov(1,1)+ delongcov(2,2)- 2*delongcov(1,2));
    pvalue_versus_shuffle(epoch,1) =  2*normcdf(abs(z), 0, 1,'upper')
        
    % one bin rate remapping VS 20ms rate fixed
    sample = [];
    sample.ratings = [log_odd.normal_zscored.rate_fixed(epoch_index{epoch});log_odd.one_bin_zscored.rate_remapping(epoch_index{epoch})];
    sample.spsizes = [length(epoch_index{epoch})/2 length(epoch_index{epoch})/2];
    [aucs, delongcov] = fastDeLong(sample);
    z = diff(aucs) / sqrt(delongcov(1,1)+ delongcov(2,2)- 2*delongcov(1,2));
    pvalue_versus_shuffle(epoch,2) =  2*normcdf(abs(z), 0, 1,'upper')
    
    % one bin rate remapping VS one bin
    sample = [];
    sample.ratings = [log_odd.one_bin_zscored.original(epoch_index{epoch}); log_odd.one_bin_zscored.rate_remapping(epoch_index{epoch})];
    sample.spsizes = [length(epoch_index{epoch})/2 length(epoch_index{epoch})/2];
    [aucs, delongcov] = fastDeLong(sample);
    z = diff(aucs) / sqrt(delongcov(1,1)+ delongcov(2,2)- 2*delongcov(1,2));
    pvalue_versus_shuffle(epoch,3) =  2*normcdf(abs(z), 0, 1,'upper')
    
    % one bin rate remapping vs rate fixed global remapping 
    sample = [];
    sample.ratings = [log_odd.normal_zscored.rate_fixed_global_remapping(epoch_index{epoch}); log_odd.one_bin_zscored.rate_remapping(epoch_index{epoch})];
    sample.spsizes = [length(epoch_index{epoch})/2 length(epoch_index{epoch})/2];
    [aucs, delongcov] = fastDeLong(sample);
    z = diff(aucs) / sqrt(delongcov(1,1)+ delongcov(2,2)- 2*delongcov(1,2));
    pvalue_versus_shuffle(epoch,4) =  2*normcdf(abs(z), 0, 1,'upper')
    
%     sample = [];
%     sample.ratings = [log_odd.normal_zscored.original(epoch_index{epoch}); log_odd.one_bin_zscored.rate_fixed(epoch_index{epoch})];
%     sample.spsizes = [length(epoch_index{epoch})/2 length(epoch_index{epoch})/2];
%     z = abs(diff(aucs)) / sqrt(delongcov(1,1)+ delongcov(2,2)- 2*delongcov(1,2));
%     pvalue_versus_shuffle(epoch,4) =  2*normcdf(abs(z), 0, 1,'upper')
end

