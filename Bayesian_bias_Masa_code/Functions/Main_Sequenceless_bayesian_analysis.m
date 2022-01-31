%% Modified Carey et al 2019 sequenceless decoding analysis
% Bayesian Decoding of Track 1 and Track 2 Replay -> Log Odd
clear all
cd D:\Rate_remapping_analysis\data;
folders = { '2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'};
[log_odd] = log_odd_bayesian_analysis(folders,1,1);

if ~isfolder({'rate_remapping_one_bin'})
    mkdir('rate_remapping_one_bin') 
end

cd rate_remapping_one_bin
save log_odd_one_bin log_odd
cd ..

clear all
cd D:\Rate_remapping_analysis\data;
folders = { '2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'};
[log_odd] = log_odd_bayesian_analysis(folders,0.02,[]);
log_odd1 = log_odd;

if ~isfolder({'rate_remapping_20_bins'})
    mkdir('rate_remapping_20_bins') 
end

cd rate_remapping_20_bins
save log_odd_20ms log_odd
cd ..

cd rate_remapping_one_bin
load log_odd_one_bin
cd ..

log_odd.normal_zscored = log_odd1.normal_zscored;
log_odd.normal = log_odd1.normal;
save log_odd log_odd
clear log_odd1

% for n = 1:5
%     cd(folders{n})  
%     load extracted_place_fields_BAYESIAN
%     for track_id = 1:2
%         place_fields_BAYESIAN.track(track_id).global_remapped_raw = [];
%         place_fields_BAYESIAN.track(track_id).random_cell = [];
%         place_fields_BAYESIAN.track(track_id).circular_shifted_raw = [];
%     end
%     save extracted_place_fields_BAYESIAN place_fields_BAYESIAN
%     cd ..
% end              

%% Wirte into table
clear all
cd D:\Rate_remapping_analysis\data;
% cd D:\Rate_remapping_analysis\backup_data\1205_best\data
% cd D:\Rate_remapping_analysis\backup_data\1105\data
load log_odd
% load ROC

folders = { '2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'};

option = {'original','rate fixed','original','rate fixed global remapping','rate remapped one bin'};
time_bin_option = [0.02,0.02,1,0.02,1];
ROC = [];


for n = 1:5
    tic
    [ROC,log_odd,AUC_combined,AUC_session,Mean_zdata,SE_zdata,pval_average,pval_session] = plot_log_odd_new(ROC,log_odd,'common good cells',time_bin_option(n),[],option{n},'Y');
    
%     [log_odd,AUC,dprime,AUC_session,dprime_session,Mean_zdata,SE_zdata,pval_average,pval_session] = plot_log_odd_new_backup(log_odd,'common good cells',time_bin_option(n),[],option{n},'Y')
%     
    toc
    close all
    label = [];
    Session(1) = {'Session'};
    %     Track_id(1) = {'Track'};
    pvalue{n}(1) = {'P value'};
    Mean{n}(1) = {'Mean Difference'};
    SE{n}(1) = {'SE'};
    behave_state(1) = {'Behavioral State'};
    Mean_AUC{n}(1) = {'Mean AUC'};
    UCI_AUC{n}(1) = {'Upper CI'};
    LCI_AUC{n}(1) = {'Lower CI'};
    Mean_dprime{n}(1) = {'Mean Dprime'};
    UCI_dprime{n}(1) = {'Upper CI'};
    LCI_dprime{n}(1) = {'Lower CI'};
    SE_AUC{n}(1) = {'SEM AUC'};
    
    c = 2;
    for epoch = 1:3
        %         for track = 1:2
        for session = 1:5
            Session(c) = {session};
            if isempty(pval_session{1}{session}{epoch})
                pvalue{n}(1,c) = {nan};              
            else
                pvalue{n}(c) = {pval_session{1}{session}{epoch}};
                Mean{n}(c) = {Mean_zdata{1}{session}(epoch)};

                SE{n}(c) = {nan};                
                    Mean_AUC{n}(c) = {mean(AUC_session{session}(epoch,:))};
                    UCI_AUC{n}(c) = {prctile(AUC_session{session}(epoch,:),97.5)};
                    LCI_AUC{n}(c) = {prctile(AUC_session{session}(epoch,:),2.5)};
                    SE_AUC{n}(c) = {std(AUC_session{session}(epoch,:))};

            end
            
            if epoch == 1

                behave_state(1,c) = {'PRE'};
            elseif epoch == 2

                behave_state(1,c) = {'RUN'};
            elseif epoch == 3

                behave_state(1,c) = {'POST'};
            end
            c = c+1;
        end

        Session(c) = {'Average'};
        Mean{n}(c) = {Mean_zdata{2}(epoch)};
        SE{n}(c) = {SE_zdata(epoch)};
        Mean_AUC{n}(c) = {mean(AUC_combined(epoch,:))};
        UCI_AUC{n}(c) = {prctile(AUC_combined(epoch,:),97.5)};
        LCI_AUC{n}(c) = {prctile(AUC_combined(epoch,:),2.5)};
        SE_AUC{n}(c) = {std(AUC_combined(epoch,:))};

        pvalue{n}(c) = {pval_average{1}{epoch}};
        behave_state(c) = {behave_state(c-1)};
        c = c +1;
    end
    %     end
end

save ROC ROC
save log_odd log_odd

cd log_odd_figure
% cd D:\Rate_remapping_analysis\data\log_odd_figure
T1 = table(behave_state',Session',Mean{1}',SE{1}',pvalue{1}',...
Mean{2}',SE{2}',pvalue{2}',Mean{3}',SE{3}',pvalue{3}',Mean{4}',SE{4}',...
pvalue{4}',Mean{5}',SE{5}',pvalue{5}')

writetable(T1,'zscore log odd difference stat.xlsx','Sheet',1,'Range','B2');


T2 = table(behave_state',Session',Mean_AUC{1}',SE_AUC{1}',UCI_AUC{1}',LCI_AUC{1}',...
Mean_AUC{2}',SE_AUC{2}',UCI_AUC{2}',LCI_AUC{2}',...
Mean_AUC{3}',SE_AUC{3}',UCI_AUC{3}',LCI_AUC{3}',...
Mean_AUC{4}',SE_AUC{4}',UCI_AUC{4}',LCI_AUC{4}',...
Mean_AUC{5}',SE_AUC{5}',UCI_AUC{5}',LCI_AUC{5}')

writetable(T2,'ROC_AUC_stat.xlsx','Sheet',1,'Range','B2');
cd ..


%% Zscore log odd distribution

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
        epoch_track_index{1}{2} = intersect(state_index,find(track_id==1));    
    elseif k == 4 % RUN Track 2
        epoch_track_index{2}{2} = intersect(state_index,find(track_id==2));        
    end
end


for epoch = 1:3
    epoch_index{epoch} = [epoch_track_index{1}{epoch} epoch_track_index{2}{epoch}];
end

data = log_odd.normal_zscored.original;

fig = figure(1)
fig.Position = [680 300 575 575];

for epoch = 1:3
    subplot(2,2,epoch)
    scatter(data(epoch_track_index{1}{epoch}),0.7+rand(1,length((epoch_track_index{1}{epoch})))*(1.3-0.7),3,'r','filled')
    hold on
    scatter(data(epoch_track_index{2}{epoch}),1.7+rand(1,length((epoch_track_index{2}{epoch})))*(2.3-1.7),3,'b','filled')
    if epoch == 1
        title('PRE')
    elseif epoch == 2
        title('RUN')
        
    elseif epoch == 3
        title('POST')
    end
    xlabel('Zscored log odd')
    yticks([1,2])
    yticklabels({'Track 1','Track 2'})
        ax = gca;
    set(ax,'LineWidth',1.5)
    ax.YAxis.TickDirection =  'out';       %repeat for XAxis
    ax.YAxis.TickLength =  [.005 1];       %repeat for XAxis
    ax.XAxis.TickDirection =  'out';       %repeat for XAxis
    ax.XAxis.TickLength =  [.005 1];       %repeat for XAxis
    ax.FontSize = 12;
end

    cd log_odd_figure
%     set(gcf,'PaperSize',[20 10]); %set the paper size to what you want
%     print(gcf,'filename','-dpdf') % then print it

    saveas(gcf,'zscore log odd distribution for original.pdf')
    cd ..


%% AUC summary

load ROC

for epoch = 1:3
    for n = 1:5
        
        if n == 1
            AUC_bootstrap{epoch}{n} = ROC(epoch).normal.original(6).AUC;
        elseif n == 2
            AUC_bootstrap{epoch}{n} = ROC(epoch).normal.rate_fixed(6).AUC;
        elseif n == 3
            AUC_bootstrap{epoch}{n} = ROC(epoch).one_bin.original(6).AUC;
        elseif n == 4
            AUC_bootstrap{epoch}{n} = ROC(epoch).normal.rate_fixed_global_remapped(6).AUC;
        elseif n == 5  
            AUC_bootstrap{epoch}{n} = ROC(epoch).one_bin.rate_remapping(6).AUC;
        end
        
    end
end

AUC = [];
AUC_LCI = [];
AUC_UCI = [];
AUC_SE = [];


for epoch = 1:3
    for n = 1:5
        AUC(epoch,n) = mean(AUC_bootstrap{epoch}{n});
%         AUC_UCI(epoch,n) = prctile(AUC_bootstrap{epoch}{n},99.375);
%         AUC_LCI(epoch,n) = prctile(AUC_bootstrap{epoch}{n},0.625);
        AUC_SE(epoch,n) = std(AUC_bootstrap{epoch}{n});
    end 
end


fig = figure(1)
fig.Position = [650 100 635 477];

color_fill = {'k','b','r'};
for epoch = 1:3
%     errorbar(X,Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1));
    p(epoch) = plot((0.5:1:4.5),AUC(epoch,:),color_fill{epoch})
    hold on
    if epoch ~= 1
        
%         scatter([1.55:1:4.55],AUC(epoch,2:end)+0.03,color_fill{epoch},'*')
    end
    hold on
    % If 95% confidence interval
%     e(epoch) = errorbar([0.5:1:4.5],AUC(epoch,:),AUC(epoch,:)-AUC_LCI(epoch,:),AUC_UCI(epoch,:)-AUC(epoch,:));
    % If SE
    e(epoch) = errorbar([0.5:1:4.5],AUC(epoch,:),AUC_SE(epoch,:),AUC_SE(epoch,:))
    e(epoch).Color = color_fill{epoch}
end

ax = gca;
set(ax,'LineWidth',1.5)
ax.YAxis.TickDirection =  'out';       %repeat for XAxis
ax.YAxis.TickLength =  [.005 1];       %repeat for XAxis
ax.XAxis.TickDirection =  'out';       %repeat for XAxis
ax.XAxis.TickLength =  [.005 1];       %repeat for XAxis
ax.FontSize = 12;


daspect([1 0.2 1])
% set(p,{'DisplayName'}, {'PRE','RUN','POST'}')
lgd = legend(p(1:3),{'PRE','RUN','POST'},'Location','northeast')
lgd.FontSize = 12;
legend boxoff

ylabel('Area under the ROC curve')
xlim([0 5])
ylim([0.3 1.1])
xticks([0.5,1.5,2.5,3.5,4.5])
xticklabels({'original','rate fixed','place/temporal removed','rate fixed place randomized','place removed rate randomized'})
xtickangle(45)
box off

cd log_odd_figure
filename = 'AUC PRE vs RUN vs POST.pdf'
set(gcf,'PaperType','A3'); %set the paper size to what you want  
print(gcf,filename,'-dpdf') % then print it
% saveas(gcf,filename)
cd ..

%% AUC Difference (VS Original)

AUC_diff_bootstrap = [];
for epoch = 1:3
    k = 1;
    for n = 2:5
        AUC_diff_bootstrap{epoch}{k} = AUC_bootstrap{epoch}{n} - AUC_bootstrap{epoch}{1};
%         AUC_diff_bootstrap{epoch}{k} = AUC_bootstrap{epoch}{1} - AUC_bootstrap{epoch}{n};
        k = k + 1;
    end
end


AUC_diff = [];
AUC_diff_LCI = [];
AUC_diff_UCI = [];
AUC_diff_SE = [];

for epoch = 1:3
    for n = 1:4
        AUC_diff(epoch,n) = mean(AUC_diff_bootstrap{epoch}{n});
%         AUC_diff_UCI(epoch,n) = prctile(AUC_diff_bootstrap{epoch}{n},99.375);
%         AUC_diff_LCI(epoch,n) = prctile(AUC_diff_bootstrap{epoch}{n},0.625);
        AUC_diff_UCI(epoch,n) = prctile(AUC_diff_bootstrap{epoch}{n},97.5);
        AUC_diff_LCI(epoch,n) = prctile(AUC_diff_bootstrap{epoch}{n},2.5);        
        
        AUC_diff_SE(epoch,n) = std(AUC_diff_bootstrap{epoch}{n});
    end 
end


fig = figure(2)
fig.Position = [650 100 635 445];
color_fill = {'k','b','r'};
b = bar((0.5:1:3.5),AUC_diff);
b(1).FaceColor = color_fill{1};
b(1).FaceAlpha = 0.50;
b(2).FaceColor = color_fill{2};
b(2).FaceAlpha = 0.50;
b(3).FaceColor = color_fill{3};
b(3).FaceAlpha = 0.50;
ax = gca;
set(ax,'LineWidth',1.5)
ax.YAxis.TickDirection =  'out';       %repeat for XAxis
ax.YAxis.TickLength =  [.005 1];       %repeat for XAxis
ax.XAxis.TickDirection =  'out';       %repeat for XAxis
ax.XAxis.TickLength =  [.005 1];       %repeat for XAxis
ax.FontSize = 12;
hold on

errorbar_xpos = [[0.28:1:3.28];[0.5:1:3.5];[0.73:1:3.73]];

for epoch = 1:3
% e(epoch) = errorbar(errorbar_xpos(epoch,:),AUC_diff(epoch,:),AUC_diff_SE(epoch,:),AUC_diff_SE(epoch,:),'.')
e(epoch) = errorbar(errorbar_xpos(epoch,:),AUC_diff(epoch,:),AUC_diff_UCI(epoch,:)-AUC_diff(epoch,:),AUC_diff(epoch,:)-AUC_diff_LCI(epoch,:),'.')

e(epoch).Color = color_fill{epoch}
end

for epoch = 2:3
s(epoch) = scatter(errorbar_xpos(epoch,1:end),AUC_diff_LCI(epoch,1:end)-0.03,color_fill{epoch},'*')
% s(3) = scatter(errorbar_xpos(3,1),AUC_diff_LCI(3,1)-0.03,color_fill{3},'*')
% % 
% s(epoch) = scatter(errorbar_xpos(epoch,2:end),AUC_diff_UCI(epoch,2:end)+0.03,color_fill{epoch},'*')
% s(3) = scatter(errorbar_xpos(3,1),AUC_diff_UCI(3,1)+0.03,color_fill{3},'*')
end


daspect([1 0.2 1])
% set(p,{'DisplayName'}, {'PRE','RUN','POST'}')
lgd = legend(b(1:3),{'PRE','RUN','POST'},'Location','southwest')
% lgd = legend(b(1:3),{'PRE','RUN','POST'},'Location','northwest')
lgd.FontSize = 12;
legend boxoff

ylabel('Area under the ROC curve difference')
xlim([0 4])
% ylim([-0.2 0.6])
ylim([-0.6 0.2])
xticks([0.5,1.5,2.5,3.5])
% xticklabels({'20ms rate fixed','one bin','20ms global remapped','one bin rate remapped'})
xticklabels({'rate fixed','place/temporal removed','rate fixed place randomized','place removed rate randomized'})
xtickangle(45)
box off


cd log_odd_figure
filename = 'AUC difference.pdf'
set(gcf,'PaperType','A3'); %set the paper size to what you want  
print(gcf,filename,'-dpdf') % then print it
% saveas(gcf,filename)
cd ..


%% AUC Difference (Vs Negative control, 20ms global remapped)

AUC_diff_bootstrap = [];
for epoch = 1:3
    k = 1;
    for n = [1 2 3 5]
%         AUC_diff_bootstrap{epoch}{k} = AUC_bootstrap{epoch}{n} - AUC_bootstrap{epoch}{1};
        AUC_diff_bootstrap{epoch}{k} = AUC_bootstrap{epoch}{n} - AUC_bootstrap{epoch}{4};
        k = k + 1;
    end
end

AUC_diff_NC = [];
AUC_diff_NC_LCI = [];
AUC_diff_NC_UCI = [];
AUC_diff_NC_SE = [];

for epoch = 1:3
    for n = 1:4
        AUC_diff_NC(epoch,n) = mean(AUC_diff_bootstrap{epoch}{n});
%         AUC_diff_NC_UCI(epoch,n) = prctile(AUC_diff_bootstrap{epoch}{n},99.375);
%         AUC_diff_NC_LCI(epoch,n) = prctile(AUC_diff_bootstrap{epoch}{n},0.625);
        AUC_diff_NC_UCI(epoch,n) = prctile(AUC_diff_bootstrap{epoch}{n},97.5);
        AUC_diff_NC_LCI(epoch,n) = prctile(AUC_diff_bootstrap{epoch}{n},2.5);
        
        AUC_diff_NC_SE(epoch,n) = std(AUC_diff_bootstrap{epoch}{n});
    end 
end


fig = figure(3)
fig.Position = [650 100 635 445];
color_fill = {'k','b','r'};
b = bar((0.5:1:3.5),AUC_diff_NC);
b(1).FaceColor = color_fill{1};
b(1).FaceAlpha = 0.50;
b(2).FaceColor = color_fill{2};
b(2).FaceAlpha = 0.50;
b(3).FaceColor = color_fill{3};
b(3).FaceAlpha = 0.50;
ax = gca;
set(ax,'LineWidth',1.5)
ax.YAxis.TickDirection =  'out';       %repeat for XAxis
ax.YAxis.TickLength =  [.005 1];       %repeat for XAxis
ax.XAxis.TickDirection =  'out';       %repeat for XAxis
ax.XAxis.TickLength =  [.005 1];       %repeat for XAxis
ax.FontSize = 12;
hold on

errorbar_xpos = [[0.28:1:3.28];[0.5:1:3.5];[0.73:1:3.73]];

for epoch = 1:3
% e(epoch) = errorbar(errorbar_xpos(epoch,:),AUC_diff(epoch,:),AUC_diff_SE(epoch,:),AUC_diff_SE(epoch,:),'.')
e(epoch) = errorbar(errorbar_xpos(epoch,:),AUC_diff_NC(epoch,:),AUC_diff_NC_UCI(epoch,:)-AUC_diff_NC(epoch,:),AUC_diff_NC(epoch,:)-AUC_diff_NC_LCI(epoch,:),'.')

e(epoch).Color = color_fill{epoch}
end

for epoch = 2:3
% s(epoch) = scatter(errorbar_xpos(epoch,2:end),AUC_diff_NC_LCI(epoch,2:end)-0.03,color_fill{epoch},'*')
% s(3) = scatter(errorbar_xpos(3,1),AUC_diff_NC_LCI(3,1)-0.03,color_fill{3},'*')
% % 
s(epoch) = scatter(errorbar_xpos(epoch,1:end-1),AUC_diff_NC_UCI(epoch,1:end-1)+0.03,color_fill{epoch},'*')
end


daspect([1 0.2 1])
% set(p,{'DisplayName'}, {'PRE','RUN','POST'}')
% lgd = legend(b(1:3),{'PRE','RUN','POST'},'Location','southwest')
lgd = legend(b(1:3),{'PRE','RUN','POST'},'Location','northeast')
lgd.FontSize = 12;
legend boxoff

ylabel('Area under the ROC curve difference')
xlim([0 4])
ylim([-0.2 0.6])
% ylim([-0.6 0.2])
xticks([0.5,1.5,2.5,3.5])
xticklabels({'Original','20ms rate fixed','one bin','one bin rate remapped'})
xticklabels({'original','rate fixed','place/temporal removed','place removed rate randomized'})

xtickangle(45)
box off


cd log_odd_figure
filename = 'AUC difference vs negative control.pdf'
set(gcf,'PaperType','A3'); %set the paper size to what you want  
print(gcf,filename,'-dpdf') % then print it
% saveas(gcf,filename)
cd ..


%% DeLong Test and confidence interval into table
load ROC
load log_odd
[pvalue_versus_original pvalue_versus_shuffle] = DeLong_test(log_odd,ROC)

T3 = table(pvalue_versus_original);
writetable(T3,'ROC_AUC_pval.xlsx','WriteVariableNames',0,'Sheet',1,'Range','G5');


T4 = table(pvalue_versus_shuffle);
writetable(T4,'ROC_AUC_pval.xlsx','WriteVariableNames',0,'Sheet',1,'Range','G12');


AUC_table = [];

for comparision = 1:2
    for epoch = 1:3
        n = 1;
        for condition = 1:4
            if comparision == 1;
                AUC_table{comparision}(epoch,n) = AUC_diff(epoch,condition);
                n = n + 1;
                AUC_table{comparision}(epoch,n) = AUC_diff_SE(epoch,condition);
                n = n + 1;
                AUC_table{comparision}(epoch,n) = AUC_diff_LCI(epoch,condition);
                n = n + 1;
                AUC_table{comparision}(epoch,n) = AUC_diff_UCI(epoch,condition);
                n = n + 1;                
            elseif comparision == 2;
                AUC_table{comparision}(epoch,n) = AUC_diff_NC(epoch,condition);
                n = n + 1;
                AUC_table{comparision}(epoch,n) = AUC_diff_NC_SE(epoch,condition);
                n = n + 1;
                AUC_table{comparision}(epoch,n) = AUC_diff_NC_LCI(epoch,condition);
                n = n + 1;
                AUC_table{comparision}(epoch,n) = AUC_diff_NC_UCI(epoch,condition);
                n = n + 1;
            end
        end
    end
end


T1 = table(AUC_table{1});
writetable(T1,'ROC_AUC_confidence.xlsx','WriteVariableNames',0,'Sheet',1,'Range','C5');


T2 = table(AUC_table{2});
writetable(T2,'ROC_AUC_confidence.xlsx','WriteVariableNames',0,'Sheet',1,'Range','C13');

%% Single event shuffled data vs actual data (for cartoon)

n = 800
data = log_odd.normal.probability_ratio(n);
shuffled_data = log_odd.normal.probability_ratio_shuffled{n};

histogram(shuffled_data,'FaceAlpha',0.5,'FaceColor','b');
hold on
plot([data data],[0 200],'r','LineWidth',5)
ylabel('Number of replay events')
xlabel('Log Odd')

