function [log_odd,AUC,dprime,AUC_session, dprime_session,Mean_zdata,SE_zdata,pval_average,pval_session] = plot_log_odd(log_odd,place_cell_type,timebin,posbin,option,figure_option,save_option)


counts = [];
% colour_line={'g','k','b','r','m'};
% colour_symbol={'go','ko','bo','ro','mo'};


states=[-1 0 1 2]; % PRE, POST and RUN


%% Replay Index for each replay (sorted by session and state)
session_index = [];
for session = 1:max(log_odd.experiment)
    for epoch = 1:length(states) % sort track id and behavioral states
        if isempty(intersect(find(log_odd.experiment == session),find(log_odd.behavioural_state == states(epoch))))
            session_index{session}{epoch} = [];
        else
            session_index{session}{epoch} = intersect(find(log_odd.experiment == session),find(log_odd.behavioural_state == states(epoch)));
        end
    end
end

track_id = log_odd.track;

index = [];% Sort replay id according to session and behavioral states
for session = 1:max(log_odd.experiment)
    for k=1:length(states) % sort track id and behavioral states
        state_index = session_index{session}{k};
        
        if isempty(state_index)
            index{session}{1}{k} = [];
            index{session}{2}{k} = [];
        else
            
            if k == 1 % PRE
                index{session}{1}{1} = intersect(state_index,find(track_id==1));
                index{session}{2}{1} = intersect(state_index,find(track_id==2));
            elseif k == 2 % POST
                index{session}{1}{3} = intersect(state_index,find(track_id==1));
                index{session}{2}{3} = intersect(state_index,find(track_id==2));
            elseif k == 3 % RUN Track 1
                index{session}{1}{2} = intersect(state_index,find(track_id==1)); % Track 1 replay on Track 1 (local replay)
            elseif k == 4 % RUN Track 2
                index{session}{2}{2} = intersect(state_index,find(track_id==2)); % Track 2 replay on Track 2 (local replay)
            end
        end
    end
end


% % Resampling so that no of T1 replay matches no of T2 replay 
% for session = 1:max(log_odd.experiment)
%     for epoch = 1:3
%         difference = abs(length(index{session}{1}{epoch}) - length(index{session}{2}{epoch})); % find the difference between two tracks
%         
%         if length(index{session}{1}{epoch}) > length(index{session}{2}{epoch})
%             if isempty(index{session}{2}{epoch})
%                 index{session}{1}{epoch} = [];
%             else
%                 index{session}{2}{epoch} = [index{session}{2}{epoch} datasample(index{session}{2}{epoch},difference)]
%             end
%             
%         elseif length(index{session}{1}{epoch}) < length(index{session}{2}{epoch})
%             if isempty(index{session}{1}{epoch}) % If track 1 is empty then track 2 also empty (happens during PRE)
%                 index{session}{2}{epoch} = [];
%             else
%                 index{session}{1}{epoch} = [index{session}{1}{epoch} datasample(index{session}{1}{epoch},difference)]
%             end
%         end
%     end
% end


%% Track ID shuffle
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
%         tempt = intersect(state_index,find(track_id==2)) % Track 2 replay on Track 1 (Preplay)
%         epoch_track_index{2}{1} = [epoch_track_index{2}{1} tempt]; % Put it into PRE              
    elseif k == 4 % RUN Track 2
%         epoch_track_index{1}{3} = intersect(state_index,find(track_id==1));
        epoch_track_index{2}{2} = intersect(state_index,find(track_id==2));        
    end
end

% Resampling so that no of T1 replay matches no of T2 replay 
% for epoch=1:3 % sort track id and behavioral states
%     difference = abs(length(epoch_track_index{1}{epoch}) - length(epoch_track_index{2}{epoch})); % find the difference between two tracks
%     combined = [epoch_track_index{1}{epoch},epoch_track_index{2}{epoch}];
%     
%     if length(epoch_track_index{1}{epoch}) > length(epoch_track_index{2}{epoch})
%         combined = [combined datasample(epoch_track_index{2}{epoch},difference)]; % combine and equate the sample size
%         balanced_index{1}{epoch} = combined(1:length(combined)/2); % divide half into one
%         balanced_index{2}{epoch} = combined(1+length(combined)/2:end); % divide another half into another one
%     elseif length(epoch_track_index{1}{epoch}) < length(epoch_track_index{2}{epoch})
%         if ~isempty(epoch_track_index{1}{epoch}) % If track 1 is empty then track 2 also empty
%             combined = [datasample(epoch_track_index{1}{epoch},difference) combined];
%             balanced_index{1}{epoch} = combined(1:length(combined)/2);
%             balanced_index{2}{epoch} = combined(1+length(combined)/2:end);
%         end
%     end
% end


% 
%% Track ID Shuffled index

shuffled_index = []; % 1000 track Id shuffles across sessions
for shuffle = 1:1000
    for epoch=1:3 % sort track id and behavioral states
        combined = [epoch_track_index{1}{epoch},epoch_track_index{2}{epoch}];
        combined = combined(randperm(length(combined)));
        shuffled_index{shuffle}{1}{epoch} = combined(1:length(epoch_track_index{1}{epoch}));
        shuffled_index{shuffle}{2}{epoch} = combined(1+length(epoch_track_index{1}{epoch}):end);
    end
end

% Resampling so that no of T1 replay matches no of T2 replay 
% shuffled_index = []; % 1000 track Id shuffles across sessions
% for shuffle = 1:1000
%     for epoch=1:3 % sort track id and behavioral states
%         difference = abs(length(epoch_track_index{1}{epoch}) - length(epoch_track_index{2}{epoch})); % find the difference between two tracks
%         combined = [epoch_track_index{1}{epoch},epoch_track_index{2}{epoch}];
%         
%         if length(epoch_track_index{1}{epoch}) > length(epoch_track_index{2}{epoch})
%             combined = [combined datasample(epoch_track_index{2}{epoch},difference)]; % combine and equate the sample size
%             combined = combined(randperm(length(combined))); % track id shuffle
%             shuffled_index{shuffle}{1}{epoch} = combined(1:length(combined)/2); % divide half into one
%             shuffled_index{shuffle}{2}{epoch} = combined(1+length(combined)/2:end); % divide another half into another one
%         elseif length(epoch_track_index{1}{epoch}) < length(epoch_track_index{2}{epoch})
%             if ~isempty(epoch_track_index{1}{epoch}) % If track 1 is empty then track 2 also empty (happens during PRE)
%                 combined = [combined datasample(epoch_track_index{1}{epoch},difference)];
%                 combined = combined(randperm(length(combined)));
%                 shuffled_index{shuffle}{1}{epoch} = combined(1:length(combined)/2);
%                 shuffled_index{shuffle}{2}{epoch} = combined(1+length(combined)/2:end);
%             end
%         end
%     end
% end
% 

% Load data
if strcmp(option,'original')
    if timebin == 1
        
        data = log_odd.one_bin_zscored.original;
%         data(isnan(data)) = 0;
        
    elseif timebin == 0.02
        
        data = log_odd.normal_zscored.original;
    end
    
elseif strcmp(option,'rate fixed')
    if timebin == 1     
        data = log_odd.one_bin_zscored.rate_fixed;
%         data(isnan(data)) = 0;
        
    elseif timebin == 0.02        
        data = log_odd.normal_zscored.rate_fixed;
    end
    
% elseif strcmp(option,'rate fixed median spike count')
%     
%     data = log_odd.one_bin_zscored.rate_fixed_median;
elseif strcmp(option,'normal circular shift')
    data = log_odd.normal_zscored.rate_fixed_circular_shift;
    
elseif strcmp(option,'rate fixed circular shift')
    data = log_odd.normal_zscored.normal_circular_shift;
elseif strcmp(option,'rate fixed shuffled spike train')
    data = log_odd.normal_zscored.rate_fixed_shuffled_spike_train;
elseif strcmp(option,'rate fixed global remapping')
    data = log_odd.normal_zscored.rate_fixed_global_remapping;
elseif strcmp(option,'rate remapped one bin') 
    data = log_odd.one_bin_zscored.rate_remapping;
end


for epoch = 1:3 % PRE, POST, RUN
    
    for nshuffle = 1:1000
        shuffle_average{epoch}(nshuffle) = mean(data(shuffled_index{nshuffle}{1}{epoch})) - mean(data(shuffled_index{nshuffle}{2}{epoch}));
    end
    

    for session= 1:max(log_odd.experiment)
        if isempty(data(index{session}{1}{epoch})) % If Preplay is missing T1 replay (ignore it)
            zdata{session}{1}{epoch} = [];
            zdata{session}{2}{epoch} = [];
        else
            
            zdata{session}{1}{epoch} = data(index{session}{1}{epoch});
            zdata{session}{2}{epoch} = data(index{session}{2}{epoch});
            
            % Difference between mean T1 log odd - T2 log odd
            Mean_zdata{1}{session}(1,epoch) = mean(zdata{session}{1}{epoch}) - mean(zdata{session}{2}{epoch});
            % I am not sure how to calculate the SE for the difference between
            % two distribution?
            
            for nshuffle = 1:1000
                tempt = [data(index{session}{1}{epoch}) data(index{session}{2}{epoch})];
                % shuffle replay event track label within the session and
                % within the behavioral epoch
                tempt = tempt(randperm(length(tempt)));
                % Calculate mean T1 (but shuffled) log odd - T2 (but shuffled) log odd
                shuffle_session{session}{epoch}(nshuffle) = mean(tempt(1:length(index{session}{1}{epoch}))) - mean(tempt(1+length(index{session}{1}{epoch}):end));
            end
            
            % Get a zscore and p value based on the distribution of the
            % event track label shuffle for each session.
            tempt = [];
            tempt = (Mean_zdata{1}{session}(1,epoch) - mean(shuffle_session{session}{epoch}))/std(shuffle_session{session}{epoch});
            
            pval_session{1}{session}{epoch} = normcdf(tempt, 0, 1,'upper'); % One tailed test
%             pval_session{2}{session}{epoch} = 2*normcdf(abs(tempt), 0, 1,'upper'); % Two tailed test (not using it)

        end       
    end
    
    tempt1 = [];
    tempt2 = [];
    for session = 1:5
        tempt1(session) = Mean_zdata{1}{session}(1,epoch);
    end
    
    Mean_zdata{2}(epoch) = mean(tempt1);
    SE_zdata(epoch) = std(tempt1)/sqrt(length(tempt1));
    
    % Get a zscore and p value based on the comparision between the real mean value and the distribution of the
    % event track label shuffle for each session.
    tempt = [];
    tempt = (Mean_zdata{2}(epoch) - mean(shuffle_average{epoch}))/std(shuffle_average{epoch});
    
    pval_average{1}{epoch} = normcdf(tempt, 0, 1,'upper'); % one tailed
%     pval_average{2}{epoch} = 2*normcdf(abs(tempt), 0, 1,'upper');
    %      pval_average{2}{epoch} = ranksum(data(balanced_index{1}{epoch})-data(balanced_index{2}{epoch}),shuffle_average{epoch});
end


if strcmp(figure_option,'Y')
    state_label = {'PRE T1-T2 log odd difference','RUN T1-T2 log odd difference','POST T1-T2 log odd difference'};
    colour_line= {'k','b','r','g'};
%     colour_line2 = {'k--','b--','r--','g--'};
    colour_symbol={'ko','bo','ro','go'};
%     colour_symbol2= {'kx','bx','rx','gx'};
    colour_fill = {[0 0 0,0.3],[0, 0, 1, 0.3],[1,0,0,0.3]};
    
    nfig = 1;
    figure(nfig)
    for epoch = 1:3
        tempt1 = [];
        tempt2 = [];
        for session= 1:5
            if ~isempty(zdata{session}{1}{epoch})
                %             scatter(mean(zdata{session}{1}{epoch}),epoch,colour_symbol{epoch});
                scatter(epoch,Mean_zdata{1}{session}(1,epoch),20,colour_symbol{epoch});
                %             tempt1 = [tempt1 mean(zdata{session}{1}{epoch})];
                hold on
                %             scatter(mean(zdata{session}{2}{epoch}),epoch,colour_symbol2{epoch});
                %             tempt2 = [tempt2 mean(zdata{session}{2}{epoch})];
                %         boxplot(epoch,zdata{session}{2}{epoch})
            end
        end
        tempt1 = Mean_zdata{2}(1,epoch);
        %     tempt2 = mean(data(epoch_track_index{2}{epoch}));
        
        h1(epoch) = plot([epoch*ones-0.25 epoch+0.25],[mean(tempt1) mean(tempt1)],colour_line{epoch},'DisplayName',state_label{epoch});
        %     h2(epoch) = plot([mean(tempt2) mean(tempt2)],[epoch*ones-0.25 epoch+0.25],colour_line2{epoch},'DisplayName',state_label2{epoch});
        %     rectangle('Position',[mean(shuffle_average{epoch})-std(shuffle_average{epoch}),epoch-0.1,std(shuffle_average{epoch})*2,0.2],'EdgeColor',colour_line{epoch});
        rectangle('Position',[epoch-0.1,mean(shuffle_average{epoch})-std(shuffle_average{epoch}),0.2,std(shuffle_average{epoch})*2],...
            'EdgeColor',colour_fill{epoch},'FaceColor', colour_fill{epoch});
        %     rectangle('Position',[mean(shuffle_average{2}{epoch})-std(shuffle_average{2}{epoch}),epoch-0.1,std(shuffle_average{2}{epoch})*2,0.2],'EdgeColor',colour_line{epoch},'LineStyle',':');
        if pval_average{1}{epoch} < 0.05
%             scatter([2.05:1:5.05],AUC(epoch,2:end)+0.03,color_fill{epoch},'*')
            scatter([epoch+0.2],[mean(tempt1)+0.1],colour_line{epoch},'*');
        end
    end
    
    
    daspect([1 1 1])
%     legend([h1],'location', 'northeast');
    % legend([h1,h2],'location', 'northeast');
    xticks([1,2,3])
    xticklabels({'PRE','RUN','POST'}) 
    ylim( [-1 4])
    xlim([0.5 3.5])
    
    ylabel('zscored log odd T1 - T2 difference')
    box off
    
    if ~isfolder({'log_odd_figure'})
        mkdir('log_odd_figure')
    end
    
    cd log_odd_figure
    filename = sprintf('log odd difference %.2i %s.pdf',timebin,option)
    saveas(gcf,filename)
    cd ..
%     
%     if ~isempty(posbin)
%         if strcmp(option,'original')
%             title({sprintf('Log Odd T1-T2 Difference - %s (%i time bin, 1 positon bins per track)',place_cell_type,timebin)})
%         else
%             title({sprintf('Log Odd T1-T2 Difference (%s) - %s (%i time bin, 1 positon bins per track)',option,place_cell_type,timebin)})
%         end
%     else
%         if strcmp(option,'original')
%             title({sprintf('Log Odd T1-T2 Difference - %s (%i time bin, 20 positon bins per track)',place_cell_type,timebin)})
%         else
%             title({sprintf('Log Odd T1-T2 Difference (%s) - %s (%i time bin, 20 positon bins per track)',option,place_cell_type,timebin)})
%         end
%     end
end


%% ROC curve
pairs = 1;

% Resampling so that no of T1 replay matches no of T2 replay 
for session = 1:max(log_odd.experiment)
    for epoch = 1:3
        difference = abs(length(index{session}{1}{epoch}) - length(index{session}{2}{epoch})); % find the difference between two tracks
        
        if length(index{session}{1}{epoch}) > length(index{session}{2}{epoch})
            if isempty(index{session}{2}{epoch})
                index{session}{1}{epoch} = [];
            else
                index{session}{2}{epoch} = [index{session}{2}{epoch} datasample(index{session}{2}{epoch},difference)]
            end
            
        elseif length(index{session}{1}{epoch}) < length(index{session}{2}{epoch})
            if isempty(index{session}{1}{epoch}) % If track 1 is empty then track 2 also empty (happens during PRE)
                index{session}{2}{epoch} = [];
            else
                index{session}{1}{epoch} = [index{session}{1}{epoch} datasample(index{session}{1}{epoch},difference)]
            end
        end
    end
end


%% If per session
for session = 1:5
    for epoch = 1:3
        session_epoch_index{session}{epoch} = [index{session}{1}{epoch} index{session}{2}{epoch}];
    end
end

index = [];
index = session_epoch_index;

track_id = zeros(1,length(log_odd.track));
track_id(find(log_odd.track == 1)) = 1;

if strcmp(figure_option,'Y')
    nfig = nfig + 1;
    figure(nfig)
    Legend = {'PRE','RUN','POST'};
end
AUC_session = [];
dprime_session = [];
FalsePR = [];
TruePR = [];
% assume track 1 is YES/1 and Track 2 is NO/0

for session = 1:5
    subplot(2,3,session)
    for epoch = 2:3
        CV_index = crossvalind('kfold',data(index{session}{epoch}),5); % divide into 5 data set indices
        
        for c = 1:max(CV_index)
            test_data = data(index{session}{epoch}(find(CV_index == c)));% Find zscored log odd 1/5 for testing
            train_data = data(index{session}{epoch}(find(CV_index ~= c)));% Find zscored log odd 4/5 for training
            test_id = track_id(index{session}{epoch}(find(CV_index == c)));% Find Track ID 1/5 for testing
            train_id = track_id(index{session}{epoch}(find(CV_index ~= c)));% Find Track ID 4/5 for training
            
            mdl = fitglm(train_data,train_id,'Distribution','binomial','Link','logit');
            scores = predict(mdl,test_data');
            
            [X,Y,T,A,OPTROCPT] = perfcurve(track_id(index{session}{epoch}(find(CV_index == c))),scores,1);% Get True and false positive rate and area under the curve
            
            [tempt,tempt_bin] = discretize(X,-0.25:0.05:1.25); % False Positive rate is discretised making X axis comparable acorss five curves
            FalsePR{session}{epoch}{c}(:) = tempt_bin(tempt);
            TruePR{session}{epoch}{c}(:) = Y;
            AUC_session{session}(epoch,c) = A;
%             if A == 1
%                 A = 0.999;
%             end
            dprime_session{session}(epoch,c) = sqrt(2)*norminv(A);

        end
        
        % calculate mean and SEM True positive rate at each X point (false positive rate)
        % (across 5 curves) this is for plotting purposes
        [sorted_FalsePR sort_index]= sort(cell2mat(FalsePR{session}{epoch}));
        tempt = cell2mat(TruePR{session}{epoch});
        sorted_TruePR = tempt(sort_index);
        mean_TruePR = splitapply(@mean,sorted_TruePR, findgroups(sorted_FalsePR));
        number_per_X = histcounts(findgroups(sorted_FalsePR),1:1:(max(findgroups(sorted_FalsePR))+1));
        SE_TruePR = splitapply(@std,sorted_TruePR, findgroups(sorted_FalsePR)) ./ sqrt(length(number_per_X));
        
        if strcmp(figure_option,'Y')
            plot(unique(sorted_FalsePR),mean_TruePR,colour_line{epoch},'LineWidth',0.5);
            hold on
            x2 = [unique(sorted_FalsePR), fliplr(unique(sorted_FalsePR))];
            inBetween = [mean_TruePR+SE_TruePR, fliplr(mean_TruePR-SE_TruePR)];
            fill(x2, inBetween, colour_line{epoch},'FaceAlpha','0.25','LineStyle','none');
            xlabel('False Positive Rate')
            ylabel('True positive rate')
            title(sprintf(['Session %i'],session))
            xlim([0 1])
            ylim([0 1])
            plot([0 1],[0 1],'k--');
            box off
        end
    end
   
end

cd log_odd_figure
filename = sprintf('ROC per session %.2i %s.pdf',timebin,option)
saveas(gcf,filename)
cd ..

%     str{epoch} = sprintf(['%s mean AUC - %.2f & mean dprime - %.2f'],Legend{epoch},mean(AUC(epoch,:)),mean(dprime(epoch,:)));
%     text(FalsePR{epoch}{c}(ceil(length(TruePR{epoch}{c})/2)),TruePR{epoch}{c}(ceil(length(TruePR{epoch}{c})/2)),str{epoch});
% if strcmp(figure_option,'Y')
%     if ~isempty(posbin)
%         if strcmp(option,'original')
%             sgtitle({sprintf('ROC for Classification by Logistic Regression - %s (%i time bin, 1 positon bins per track)',place_cell_type,timebin)})
%         else
%             sgtitle({sprintf('ROC for Classification by Logistic Regression(%s) - %s (%i time bin, 1 positon bins per track)',option,place_cell_type,timebin)})
%         end
%     else
%         if strcmp(option,'original')
%             sgtitle({sprintf('ROC for Classification by Logistic Regression - %s (%i time bin, 20 positon bins per track)',place_cell_type,timebin)})
%         else
%             sgtitle({sprintf('ROC for Classification by Logistic Regression(%s) - %s (%i time bin, 20 positon bins per track)',option,place_cell_type,timebin)})
%         end
%     end
% end


%% if session combined

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


epoch_index = [];
for epoch = 1:3
    epoch_index{epoch} = [balanced_index{1}{epoch} balanced_index{2}{epoch}];
end

% epoch_index = [];
% for k=1:length(states)
%     state_index = find(log_odd.behavioural_state==states(k));
%     
%     if states(k) == 1 % On track 1
%         epoch_index{2} = [intersect(state_index,find(log_odd.track==1))]; % Track 1 replay on Track 1 considered RUN Replay
%     elseif states(k) == 2 % On track 2
%         tempt = intersect(state_index,find(log_odd.track==2));
%          epoch_index{2} =  [epoch_index{2} tempt];
%     elseif k == 2
%         epoch_index{3} = [intersect(state_index,find(log_odd.track==1)) intersect(state_index,find(log_odd.track==2))];
%     elseif k == 1
%         epoch_index{1} = [intersect(state_index,find(log_odd.track==1)) intersect(state_index,find(log_odd.track==2))];
%     end
% end

index = epoch_index;
track_id = zeros(1,length(log_odd.track));
track_id(find(log_odd.track == 1)) = 1;


if strcmp(figure_option,'Y')
    nfig = nfig + 1;
    figure(nfig)
    Legend = {'PRE','RUN','POST'};
end

AUC = [];
dprime = [];
FalsePR = [];
TruePR = [];
% assume track 1 is YES/1 and Track 2 is NO/0


for epoch = 1:3
    CV_index = crossvalind('kfold',data(index{epoch}),5); % divide into 5 roughly equal dataset
    length(index{epoch}) - (floor(length(index{epoch})/5)*5)
    
    for c = 1:max(CV_index)
        test_data = data(index{epoch}(find(CV_index == c)));% Find zscored log odd 1/5 for testing
        train_data = data(index{epoch}(find(CV_index ~= c)));% Find zscored log odd 4/5 for training
        test_id = track_id(index{epoch}(find(CV_index == c)));% Find Track ID 1/5 for testing
        train_id = track_id(index{epoch}(find(CV_index ~= c)));% Find Track ID 4/5 for training
        
        mdl = fitglm(train_data,train_id,'Distribution','binomial','Link','logit');
        scores = predict(mdl,test_data');
        
        
        [X,Y,T,A,OPTROCPT] = perfcurve(track_id(index{epoch}(find(CV_index == c))),scores,1);% Get True and false positive rate and area under the curve
        [tempt,tempt_bin] = discretize(X,-0.25:0.05:1.25); % False Positive rate is discretised making X axis comparable acorss five curves
        %         tempt_bin = tempt_bin + 0.05;
        FalsePR{epoch}{c}(:) = tempt_bin(tempt);
        
        TruePR{epoch}{c}(:) = Y;
        AUC(epoch,c) = A;
        
        dprime(epoch,c) = sqrt(2)*norminv(A);
        
    end
    
    % calculate mean and SEM True positive rate at each false positive rate bin
    % (across 5 curves) this is for plotting purposes
    [sorted_FalsePR sort_index]= sort(cell2mat(FalsePR{epoch}));
    tempt = cell2mat(TruePR{epoch});
    sorted_TruePR = tempt(sort_index);
    mean_TruePR = splitapply(@mean,sorted_TruePR, findgroups(sorted_FalsePR));
    number_per_X = histcounts(findgroups(sorted_FalsePR),1:1:(max(findgroups(sorted_FalsePR))+1));
    SE_TruePR = splitapply(@std,sorted_TruePR, findgroups(sorted_FalsePR)) ./ sqrt(length(number_per_X));
    
    
    % Save ROC curve
    if timebin == 0.02
        if strcmp(option,'original')
            log_odd.ROC(epoch).normal.original.sorted_FalsePR = unique(sorted_FalsePR);
            log_odd.ROC(epoch).normal.original.SE_TruePR = SE_TruePR;
            log_odd.ROC(epoch).normal.original.mean_TruePR = mean_TruePR;
        elseif strcmp(option,'rate fixed')
            log_odd.ROC(epoch).normal.rate_fixed.sorted_FalsePR = unique(sorted_FalsePR);
            log_odd.ROC(epoch).normal.rate_fixed.SE_TruePR = SE_TruePR;
            log_odd.ROC(epoch).normal.rate_fixed.mean_TruePR = mean_TruePR;
        elseif strcmp(option,'rate fixed global remapping')
            log_odd.ROC(epoch).normal.rate_fixed_global_remapped.sorted_FalsePR = unique(sorted_FalsePR);
            log_odd.ROC(epoch).normal.rate_fixed_global_remapped.SE_TruePR = SE_TruePR;
            log_odd.ROC(epoch).normal.rate_fixed_global_remapped.mean_TruePR = mean_TruePR;
        end
        
    elseif timebin == 1
        if strcmp(option,'original')
            log_odd.ROC(epoch).one_bin.original.sorted_FalsePR = unique(sorted_FalsePR);
            log_odd.ROC(epoch).one_bin.original.SE_TruePR = SE_TruePR;
            log_odd.ROC(epoch).one_bin.original.mean_TruePR = mean_TruePR;
        elseif strcmp(option,'rate fixed')
            log_odd.ROC(epoch).one_bin.rate_fixed.sorted_FalsePR = linspace(0,1,21);
            log_odd.ROC(epoch).one_bin.rate_fixed.SE_TruePR = linspace(SE_TruePR(1),SE_TruePR(end),21);
            log_odd.ROC(epoch).one_bin.rate_fixed.mean_TruePR = linspace(mean_TruePR(1),mean_TruePR(end),21);
        end
    end
    
    
    if strcmp(figure_option,'Y')
        p(epoch)= plot(unique(sorted_FalsePR),mean_TruePR,colour_line{epoch});
        hold on
        x2 = [unique(sorted_FalsePR), fliplr(unique(sorted_FalsePR))];
        inBetween = [mean_TruePR+SE_TruePR, fliplr(mean_TruePR-SE_TruePR)];
        fill(x2, inBetween, colour_line{epoch},'FaceAlpha','0.25','LineStyle','none');
        xlabel('False Positive Rate')
        ylabel('True positive rate')
        
%         str{epoch} = sprintf(['mean AUC - %.3f\nmean dprime - %.3f'],mean(AUC(epoch,:)),mean(dprime(epoch,:)));
% %         str{epoch} = sprintf(['%s\nmean AUC - %.3f\nmean dprime - %.3f'],Legend{epoch},mean(AUC(epoch,:)),mean(dprime(epoch,:)));
%         if epoch == 1
%             text(FalsePR{epoch}{c}(ceil(length(TruePR{epoch}{c})/2)),TruePR{epoch}{c}(ceil(length(TruePR{epoch}{c})/2)),str{epoch},'Color',colour_line{epoch});
%         elseif epoch == 2
%             text(FalsePR{epoch}{c}(ceil(length(TruePR{epoch}{c})/2)),[TruePR{epoch}{c}(ceil(length(TruePR{epoch}{c})/2)) - 0.05],str{epoch},'Color',colour_line{epoch});
%         elseif epoch == 3
%             text(FalsePR{epoch}{c}(ceil(length(TruePR{epoch}{c})/2)),[TruePR{epoch}{c}(ceil(length(TruePR{epoch}{c})/2))- 0.05],str{epoch},'Color',colour_line{epoch});
%         end
%         
       
        xlim([0 1])
        ylim([0 1])
        plot([0 1],[0 1],'k--');
        box off
    end
     
end
legend(p(1:3),{'PRE','RUN','POST'})

cd log_odd_figure
filename = sprintf('ROC session combined %.2i %s.pdf',timebin,option)
saveas(gcf,filename)
cd ..

% 
% if strcmp(figure_option,'Y')
%     if ~isempty(posbin)
%         if strcmp(option,'original')
%             sgtitle({sprintf('ROC for Classification by Logistic Regression - %s (%i time bin, 1 positon bins per track)',place_cell_type,timebin)})
%         else
%             sgtitle({sprintf('ROC for Classification by Logistic Regression(%s) - %s (%i time bin, 1 positon bins per track)',option,place_cell_type,timebin)})
%         end
%     else
%         if strcmp(option,'original')
%             sgtitle({sprintf('ROC for Classification by Logistic Regression - %s (%i time bin, 20 positon bins per track)',place_cell_type,timebin)})
%         else
%             sgtitle({sprintf('ROC for Classification by Logistic Regression(%s) - %s (%i time bin, 20 positon bins per track)',option,place_cell_type,timebin)})
%         end
%     end
%     nfig = nfig + 1;
% end
% 



end
