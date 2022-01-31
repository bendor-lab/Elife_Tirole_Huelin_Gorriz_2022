%% use best sessions
load('folders_to_process_REPLAY.mat');
fd= folders([8 13 14 16 17]);
remapping= rate_remapping_TRACK_PAIRS(fd,'wcorr',0);
plot_rate_remapping_NEW('use_mat',remapping,...
                        'subset',{'stable cells laps'},'epoch',{'PRE','POST'},...
                         'x_label','peak in-field rate difference (Hz)',...
                         'y_label','peak replay rate difference (Hz)');
                     
%% cell number resampling
plot_rate_remapping_NEW('subset',{'stable cells laps'},'epoch',{'PRE','POST'},...
                         'x_label','peak in-field rate difference (Hz)',...
                         'y_label','peak replay rate difference (Hz)');

num_to_resample= 124; % from stable_cells in original data
load('U:\Margot\Reward_Experiment\Tables\subsets_of_cells.mat');
subset_of_cells.subset{2}= 'rand_resampling';
for i=1:1000
    for track_pair=1:length(subset_of_cells.cell_IDs{1})
        avail_cells= subset_of_cells.cell_IDs{1}{track_pair};
        rnd_idx= randi(length(avail_cells),num_to_resample,1);
        subset_of_cells.cell_IDs{2}{track_pair}= avail_cells(rnd_idx);
    end
    save('./Tables/subsets_of_cells.mat','subset_of_cells');
    % already a subset of stable cells
    [pval,F,~]= plot_rate_remapping_NEW('subset',{'rand_resampling'},'epoch',{'PRE','POST'});
    close all
    
    resampling(i).pval_PRE= [pval(1,:,1) pval(1,:,2) pval(1,:,3)];
    resampling(i).pval_POST= [pval(2,:,1) pval(2,:,2) pval(2,:,3)];
end

save('.\Tables\reg_resampled.mat','resampling');

all_pval_PRE= [resampling.pval_PRE];
all_pval_POST= [resampling.pval_POST];
step1= 0.005; step2= 0.01;
edg1= 0:step1:0.05; edg2= 0.06:step2:0.5;
edg= [edg1 edg2];
ctr= [edg1(1:end-1)+step1./2 0.055 edg2(1:end-1)+step2./2];
figure('Color','w');
subplot(2,2,1);
n= histcounts(all_pval_PRE,edg,'Normalization','countdensity');
bar(ctr,n);hold on;
plot([0.05 0.05],ylim,'r');
xlabel('pval');
ylabel('count density');
title('PRE')
subplot(2,2,2);
edg_post= [0:0.0001:0.005];
n= histcounts(all_pval_POST,edg_post,'Normalization','countdensity');
bar(edg_post(1:end-1)+0.0001./2,n);
xlabel('pval');
title('POST')
