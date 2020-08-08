function figs = dm_measurestatplot(cfg,data,stats)
% dm_measurestatplot plots the output of ft_measurestatistics or
% ft_applymeasure (if stats input is empty)
%
% Inputs:
%
% cfg: a structure with the following fields:
%      datatype: 'eeg', 'meg', or 'source' (default = 'eeg')
%      cond: names for the conditions to be plotted. Should be a cell
%         array of the same length as data (default = {'Condition
%         1','Condition 2', etc}
%      channel: a channel input that works with ft_channelselection
%         (default = 'all')
%      meas: indices of the measures you want to plot (default = 'all' -
%         makes a new figure for each measure
%      measname: a cell array indicating the names of the measures
%         (default =from measure handles)
%      lay: a fieldtrip layout if using meg data
%      plotmode: 'topo','violin', or 'combined'. Plots either the topography
%         of the measures or a violin plot of the means for each group.
%         'combined' plots the violin plots below and the topoplots above
%         (default = 'combined')
%      colormap: the colormap to use when plotting topographies (default =
%         parula)
%      plotparam: the parameter to be plotted topographically for the
%         difference in conditions - either 'dif' for the difference or
%         'effsize' for the effect size statistic (default = 'dif')
%      datasetinfo: used for source-space plotting. Should contain the
%         fields:
%
%         atlas: the atlas used to parcellate the data
%         sourcemodel: a head shape or source model that corresponds to the
%         atlas
%      savefig: save the figures as Matlab figures (default = 'no')
%      reverse: if set to 1, reverses the comparison condition from cond_1 
%         - cond_2 to cond_2 - cond_1 (default = 0) 
%
% data: a cell array of outputs structs from ft_applymeasure
%
% stats: a stats cell array from ft_measurestatistics. If no input is
%      given, only the topography of the data is plotted

%% Set up defaults
if ~cfgcheck(cfg,'datatype')
    cfg.datatype = 'eeg';
end

if ~cfgcheck(cfg,'cond')
    cfg.cond = repmat({'Condition'},1,length(data));
    for c = 1:length(cfg.cond)
        cfg.cond{c} = [cfg.cond{c} '_' num2str(c)];
    end
end

if ~cfgcheck(cfg,'channel')
    cfg.channel = 1:size(data{1}.data,2);
end

if ~cfgcheck(cfg,'meas')
    cfg.meas = 1:length(data{1}.meas);
end

cfg = setdefault(cfg,'reverse',0);

cfg = setdefault(cfg,'measname',cellfun(@func2str,data{1}.meas(cfg.meas),'UniformOutput',false));

cfg = setdefault(cfg,'colormap',parula);

if ~cfgcheck(cfg,'plotmode')
    cfg.plotmode = 'combined';
end

cfg = setdefault(cfg,'plotparam','dif');

cfg = setdefault(cfg,'savefig','no');

if nargin < 3
    stats = [];
end

if cfg.reverse
    cord = [2 1];
else
    cord = [1 2];
end

if isfield(data{1},'elec')
    chans = ft_channelselection(cfg.channel,data{1}.elec);
    chans = match_str(data{1}.chan,chans);
elseif isfield(data{1},'grad')
    chans = ft_channelselection(cfg.channel,data{1}.grad);
    chans = match_str(data{1}.chan,chans);
elseif iscell(cfg.channel)
    chans = match_str(data{1}.chan,cfg.channel);
else  %assuming index of channels
    chans = cfg.channel;
end

%% Plotting topos

if cfgcheck(cfg,'plotmode','topo')
    if ~isempty(stats)
        if cfgcheck(cfg,'datatype','eeg')
            for c = cfg.meas
                figs(c) = figure;
                clear ax
                for cc = 1:length(data)
                    subplot(1,length(data)+1,cc)
                    topoplot(mean(data{cc}.data(:,chans,c),1),data{1}.chanlocs(chans),'maplimits','maxmin');
                    title(cfg.cond{cc})
                    FixAxes(gca,20)
                    colormap(cfg.colormap)
                    ax(cc) = gca;
                end
                Normalize_Clim(ax)
                cbar = colorbar('Location','eastoutside'); cbar.Label.String = cfg.measname{find(cfg.meas==c)}; cbar.Label.FontSize = 20;
                
                
                subplot(1,length(data)+1,length(data)+1)
                if cfgcheck(cfg,'plotparam','dif')
                    plotparam = mean(data{cord(1)}.data(:,chans,c),1)-mean(data{cord(2)}.data(:,chans,c),1);
                elseif cfgcheck(cfg,'plotparam','effsize')
                    tmp = cat(1,stats{c}.effsize{:});
                    plotparam = horz(extractfield(tmp,stats{c}.cfg.effectsize));
                end
                
                if isfield(stats{c},'cluster')
                    cluster_topoplot(plotparam,...
                        data{1}.chanlocs(chans),stats{c}.p,stats{c}.cluster.mask)
                elseif isfield(stats{c},'fdr')
                    cluster_topoplot(plotparam,...
                        data{1}.chanlocs(chans),stats{c}.p,stats{c}.fdr < 0.05)
                else
                    cluster_topoplot(plotparam,...
                        data{1}.chanlocs(chans),stats{c}.p,stats{c}.p < 0.05)
                end
                title([cfg.cond{cord(1)} ' - ' cfg.cond{cord(2)}])
                FixAxes(gca,20)
                colormap(cfg.colormap)
                cbar = colorbar('Location','eastoutside'); cbar.Label.String = ['Difference in ' cfg.measname{find(cfg.meas==c)}]; cbar.Label.FontSize = 20;
                set(figs(c),'name',cfg.measname{find(cfg.meas==c)},'color','w')
            end
        elseif cfgcheck(cfg,'datatype','meg')
            for c = cfg.meas
                figs(c) = figure;
                clear ax
                for cc = 1:length(data)
                    subplot(1,length(data)+1,cc)
                    ft_topoplot_vec(cfg.lay,mean(data{cc}.data(:,chans,c),1),data{1}.chan(chans));
                    title(cfg.cond{cc})
                    FixAxes(gca,20)
                    colormap(cfg.colormap)
                    ax(cc) = gca;
                end
                Normalize_Clim(ax)
                
                if cfgcheck(cfg,'plotparam','dif')
                    plotparam = mean(data{cord(1)}.data(:,chans,c),1)-mean(data{cord(2)}.data(:,chans,c),1);
                elseif cfgcheck(cfg,'plotparam','effsize')
                    tmp = cat(1,stats{c}.effsize{:});
                    plotparam = horz(extractfield(tmp,stats{c}.cfg.effectsize));
                end
                
                subplot(1,length(data)+1,length(data)+1)
                if isfield(stats{c},'cluster')
                    ft_cluster_topoplot(cfg.lay,plotparam,...
                        data{1}.chan(chans),stats{c}.p,stats{c}.cluster.mask)
                elseif isfield(stats{c},'fdr')
                    ft_cluster_topoplot(cfg.lay,plotparam,...
                        data{1}.chan(chans),stats{c}.p,stats{c}.fdr < 0.05)
                else
                    ft_cluster_topoplot(cfg.lay,plotparam,...
                        data{1}.chan(chans),stats{c}.p,stats{c}.p < 0.05)
                end
                title([cfg.cond{cord(1)} ' - ' cfg.cond{cord(2)}])
                cbar = colorbar('Location','eastoutside'); cbar.Label.String = ['Difference in ' cfg.measname{find(cfg.meas==c)}]; cbar.Label.FontSize = 20;
                FixAxes(gca,20)
                colormap(cfg.colormap)
                set(figs(c),'name',cfg.measname{find(cfg.meas==c)},'color','w')
                
            end
            
        elseif cfgcheck(cfg,'datatype','source')
            for c = cfg.meas
                figs(c) = figure;
                clear ax
                for cc = 1:length(data)
                    subplot(1,length(data)+1,cc)
                    ft_cluster_sourceplot(mean(data{cc}.data(:,chans,c),1),cfg.datasetinfo.sourcemodel,cfg.datasetinfo.atlas,...
                        ones(size(mean(data{cc}.data(:,chans,c),1)))); % won't work! Fix later
                    title(cfg.cond{cc})
                    FixAxes(gca,20)
                    colormap(cfg.colormap)
                    ax(cc) = gca;
                end
                Normalize_Clim(ax)
                
                if cfgcheck(cfg,'plotparam','dif')
                    plotparam = mean(data{cord(1)}.data(:,chans,c),1)-mean(data{cord(2)}.data(:,chans,c),1);
                elseif cfgcheck(cfg,'plotparam','effsize')
                    tmp = cat(1,stats{c}.effsize{:});
                    plotparam = horz(extractfield(tmp,stats{c}.cfg.effectsize));
                end
                
                subplot(1,length(data)+1,length(data)+1)
                if isfield(stats{c},'cluster')
                    ft_cluster_sourceplot(plotparam,...
                        cfg.datasetinfo.sourcemodel,cfg.datasetinfo.atlas,stats{c}.cluster.mask)
                elseif isfield(stats{c},'fdr')
                    ft_cluster_topoplot(plotparam,...
                        cfg.datasetinfo.sourcemodel,cfg.datasetinfo.atlas,stats{c}.fdr < 0.05)
                else
                    ft_cluster_topoplot(cfg.lay,plotparam,...
                        cfg.datasetinfo.sourcemodel,cfg.datasetinfo.atlas,stats{c}.p < 0.05)
                end
                title([cfg.cond{cord(1)} ' - ' cfg.cond{cord(2)}])
                cbar = colorbar('Location','eastoutside'); cbar.Label.String = ['Difference in ' cfg.measname{find(cfg.meas==c)}]; cbar.Label.FontSize = 20;
                %FixAxes(gca,20)
                colormap(cfg.colormap)
                set(figs(c),'name',cfg.measname{find(cfg.meas==c)},'color','w')
                
            end
            
        end
    else
        if cfgcheck(cfg,'datatype','eeg')
            for c = cfg.meas
                figs(c) = figure;
                for cc = 1:length(data)
                    subplot(1,length(data),cc)
                    topoplot(mean(data{cc}.data(:,chans,c),1),data{1}.chanlocs(chans),'maplimits','maxmin');
                    title(cfg.cond{cc})
                    FixAxes(gca,20)
                    colormap(cfg.colormap)
                end
                Normalize_Clim(gcf)
                cbar = colorbar('Location','eastoutside'); cbar.Label.String = cfg.measname{find(cfg.meas==c)}; cbar.Label.FontSize = 20;
                set(figs(c),'name',cfg.measname{find(cfg.meas==c)},'color','w')
                
            end
            
        elseif cfgcheck(cfg,'datatype','meg')
            for c = cfg.meas
                figs(c) = figure;
                for cc = 1:length(data)
                    subplot(1,length(data),cc)
                    ft_topoplot_vec(cfg.lay,mean(data{cc}.data(:,chans,c),1),data{1}.chan(chans));
                    title(cfg.cond{cc})
                    FixAxes(gca,20)
                    colormap(cfg.colormap)
                end
                Normalize_Clim(gcf)
                cbar = colorbar('Location','eastoutside'); cbar.Label.String = cfg.measname{find(cfg.meas==c)}; cbar.Label.FontSize = 20;
                set(figs(c),'name',cfg.measname{find(cfg.meas==c)},'color','w')
                
            end
        elseif cfgcheck(cfg,'datatype','source')
            for c = cfg.meas
                figs(c) = figure;
                for cc = 1:length(data)
                    subplot(1,length(data)+1,cc)
                    ft_cluster_sourceplot(mean(data{cc}.data(:,chans,c),1),cfg.datasetinfo.sourcemodel,cfg.datasetinfo.atlas,...
                        ones(size(mean(data{cc}.data(:,chans,c),1)))); %won't work! Fix later
                    title(cfg.cond{cc})
                    FixAxes(gca,20)
                    colormap(cfg.colormap)
                    Normalize_Clim(gcf,1)
                end
                cbar = colorbar('Location','eastoutside'); cbar.Label.String = cfg.measname{find(cfg.meas==c)}; cbar.Label.FontSize = 20;
                set(figs(c),'name',cfg.measname{find(cfg.meas==c)},'color','w')
                
            end
        end
    end
elseif cfgcheck(cfg,'plotmode','violin')
    %% Plotting violins
    for i = cfg.meas
        figs(i) = figure;
        for c = 1:length(data)
            datastruct.(cfg.cond{c}) = mean(data{c}.data(:,chans,i),2);
        end
        ylabel(cfg.measname{find(cfg.meas==i)})
        violinplot(datastruct)
        FixAxes(gca,20)
        set(figs(i),'name',cfg.measname{find(cfg.meas==i)},'color','w');
    end
elseif cfgcheck(cfg,'plotmode','combined')
    for i = cfg.meas
        tmpcfg = cfg; tmpcfg.plotmode = 'topo'; tmpcfg.meas = i; tmpcfg.savefig = 'no'; tmpcfg.channel = chans;
        topofig = dm_measurestatplot(tmpcfg,data,stats);
        topofig = topofig(i);
        
        tmpcfg.plotmode = 'violin'; tmpcfg.meas = i; tmpcfg.savefig = 'no'; tmpcfg.channel = chans;
        violinfig = dm_measurestatplot(tmpcfg,data,stats);
        violinfig = violinfig(i);
        
        figs(i) = figure;
        p = panel('no-manage-font');
        p.pack('v',{40 60})
        nplots = length(data)+(~isempty(stats));
        p(1).pack('h',[repmat({96/nplots},1,nplots) {4}]);
        figaxes = findobj('Parent',topofig,'Type','axes');
        for c = 1:nplots
            p(1,c).select(figaxes(nplots-c+1));
        end
        %cbar = findobj('Parent',topofig,'Type','colorbar');
        %p(1,length(data)+1).select(cbar);
        colormap(cfg.colormap)
        
        
        figaxes = findobj('Parent',violinfig,'Type','axes');
        p(2).select(figaxes)
        
        p.margintop = 10;
        p.marginleft = 25;
        p(1).marginbottom = 5;
        
        %AddFigureLabel(p(1,1).axis,'A','yes')
        %AddFigureLabel(p(2).axis,'B')
        
        set(figs(i),'name',cfg.measname{find(cfg.meas==i)},'color','w','units','normalized','position',[0.15 0.3 0.7 0.7])
        
        close(topofig)
        close(violinfig)
    end
end

if cfgcheck(cfg,'savefig','yes')
    for c = 1:length(figs)
        savefig(figs(c),['Fig ' num2str(c) '-' cfg.measname{c} '.fig']);
    end
end