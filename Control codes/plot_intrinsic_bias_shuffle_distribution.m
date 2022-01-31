function plot_intrinsic_bias_shuffle_distribution(file_location)
% options: 'rate_bias' for the track peak rate shuffle
 %                   'replay_rate_bias' for the replay peak rate shuffle

% Load real correlation
[~,F_real] = plot_rate_remapping_NEW('x_var',{'place_field_diff'},'y_var',{'mean_max_FR_replay_diff'},...
                                                                                'epochs',{'PRE','POST'},'subset','stable cells laps');
close(gcf)

% Load shuffle distribution & correlation
if strcmp(file_location,'rate_bias')
    load([pwd '\CONTROLS\rate_remapped\intrinsic_rate_bias_shuffle_dist.mat'])
    fig_name=  'Track peak rate shuffle - Intrinsic bias';
    bias_shuffle_dist= intrinsic_rate_bias_shuffle_dist;
    [pval,F] = plot_rate_remapping_NEW('x_var',{'place_field_diff'},'y_var',{'mean_max_FR_replay_diff'},'epochs',{'PRE','POST'},...
                                                            'control','rate_intrinsic_bias_control','subset','stable cells laps',...
                                                            'x_label','Peak Rate Change (Hz)','y_label','Replay Rate Change (Hz)');
elseif strcmp(file_location,'replay_rate_bias')
    load([pwd '\CONTROLS\replay_rate_shuffle\intrinsic_replay_rate_bias_shuffle_dist.mat'])
    fig_name=  'Replay peak rate shuffle - Intrinsic bias';
    bias_shuffle_dist= intrinsic_replay_rate_bias_shuffle_dist;
     [pval,F] = plot_rate_remapping_NEW('x_var',{'place_field_diff'},'y_var',{'mean_max_FR_replay_diff'},'epochs',{'PRE','POST'},...
                                                            'control','replay_rate_shuffle_control','subset','stable cells laps',...
                                                            'x_label','Peak Rate Change (Hz)','y_label','Replay Rate Change (Hz)');
end


f1=gcf;
ax3 = flipud(f1.Children);

f2 = figure('units','normalized','Color','w');
f2.Name = fig_name;

for thisax = 1 : 2
    ax(thisax) = subplot(2,1,thisax);
    hold on

    h = histogram(ax(thisax),log(bias_shuffle_dist.Fstat(:,thisax)),'BinWidth',0.3,'EdgeColor','k','FaceColor',[0 0 0],'FaceAlpha',1,'EdgeAlpha',1);
    lm = prctile(log(bias_shuffle_dist.Fstat(:,thisax)),95);
    plot(ax(thisax),[lm lm],[min(ylim) max(ylim)],'LineWidth',2,'LineStyle',':','Color',[0.6 0.6 0.6])
    plot(ax(thisax),[log(F_real(thisax)) log(F_real(thisax))],[min(ylim) max(ylim)],'LineWidth',2.5,'Color',[0.5 0 0])
    xlabel(ax(thisax),'log(F-statistic)')
    ylabel(ax(thisax),'Shuffle count')
    xtickformat('%i');
    ytickformat('%i');
    yticks(ax,[25 50 75 100]);
    
    if thisax ==1    
        ax2(thisax) = axes('Position',[.22 .77 .25 .15]);        
        chi = get(ax3(thisax),'children');
        copyobj(chi,ax2(thisax));        
        annotation('textbox',[.39 .83 .1 .1],'String',{['F = ' num2str(round(log(F(thisax)),2))]; ['p = ' num2str(round(pval(thisax),2))]},'EdgeColor','none','FontSize',10)
        ax2(thisax).YLim = ax3(thisax).YLim; 
        ax2(thisax).YTick = ax3(thisax).YTick;
        ax2(thisax).YTickLabel = ax3(thisax).YTickLabel;
        ax2(thisax).XLim = ax3(thisax).XLim; 
        ax2(thisax).XTick = ax3(thisax).XTick; 
        ax2(thisax).XTickLabel = ax3(thisax).XTickLabel;
        ylh = ylabel(ax2(thisax),ax3(thisax).YAxis.Label.String); 
        xlabel(ax2(thisax),ax3(thisax).XAxis.Label.String);
        xtickformat('%i');
    
    else
        %break_axis('handle',ax(thisax),'axis','x','length',30,'start_position',15,'end_position',20) %axis break
        
        ax2(thisax) = axes('Position',[.22 .3 .25 .15]);
        chi = get(ax3(thisax),'children');
        copyobj(chi,ax2(thisax));
        annotation('textbox',[.39 .36 .1 .1],'String',{['F = ' num2str(round(log(F(thisax)),2))]; ['p = ' num2str(round(pval(thisax),2))]},'EdgeColor','none','FontSize',10)
        ax2(thisax).YLim = ax3(thisax).YLim; 
        ax2(thisax).YTick = ax3(thisax).YTick;
        ax2(thisax).YTickLabel = ax3(thisax).YTickLabel;
        ax2(thisax).XLim = ax3(thisax).XLim; 
        ax2(thisax).XTick = ax3(thisax).XTick;
        ax2(thisax).XTickLabel = ax3(thisax).XTickLabel;
        ylh = ylabel(ax2(thisax),ax3(thisax).YAxis.Label.String); 
        xlabel(ax2(thisax),ax3(thisax).XAxis.Label.String);
        xtickformat('%i');
    end


end

run_format_settings(f2,'match_ax')
end

