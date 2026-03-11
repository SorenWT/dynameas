function tbl = dm_measurestatsummary(cfg,stats,data)
% ft_statsummary outputs a summary table of the results obtained from
% ft_measurestatistics.
%
% Input arguments:
%
% cfg: a config structure with the following fields:
%      cond: the names of the conditions (default =
%         {'Condition1','Condition2'}
%      thresh: the significance threshold for including a result in the
%         table (default = 0.1)
%      test: whether to use the mean or median to determine difference
%         (default = based on cfg.test in the stats structure)
% stats: the output from dm_measurestatistics
%
% data: a cell array of outputs from dm_applymeasure - these should be in
%      the SAME ORDER as input to dm_measurestatistics. Optional input: if
%      not input, data will be taken from stats
%
% Outputs:
%
% tbl: a summary table including the measure name, p value, labels included
%      in the cluster, mean values across electrodes, mean values across
%      the cluster, difference within the cluster, and whether the measure
%      was higher in cfg.cond{1} or cfg.cond{2}


if isempty(cfg)
    cfg = struct;
end

if nargin < 3 || isempty(data)
    data = stats{1}.data;  
end

cfg = setdefault(cfg,'thresh',0.2);
cfg = setdefault(cfg,'cond',{'Condition1','Condition2'});
if ~isfield(cfg,'test')
    switch stats{1}.cfg.test
        case {'ranksum','signrank'}
            cfg.test = @median;
        case {'ttest','ttest2'}
            cfg.test = @mean;
    end
end

count = 1;
for c = 1:length(stats)
    if isfield(stats{c}.cluster,'posclusters')
        for cc = 1:length(stats{c}.cluster.posclusters)
            if stats{c}.cluster.posclusters(cc).prob < cfg.thresh
                try
                    measure{count} = func2str(data{1}.meas{c});
                catch
                    measure{count} = num2str(c);
                end
                pvalue(count) = stats{c}.cluster.posclusters(cc).prob;
                tmp = stats{c}.cluster.label(find(stats{c}.cluster.posclusterslabelmat == cc));
                tmp = cellcat(', ',tmp,'',1);
                labels{count} = [tmp{:}];
                labels{count}(end-1:end) = [];
                if ~isempty(data)
                    allmean1(count) = cfg.test(mean(data{1}.data(:,:,c),2));
                    allmean2(count) = cfg.test(mean(data{2}.data(:,:,c),2));
                    
                    clustmean1(count) = cfg.test(mean(data{1}.data(:,find(stats{c}.cluster.posclusterslabelmat == cc),c),2));
                    clustmean2(count) = cfg.test(mean(data{2}.data(:,find(stats{c}.cluster.posclusterslabelmat == cc),c),2));
                    
                    diff(count) = clustmean1(count)-clustmean2(count);
                    if diff(count) > 0
                        dirn{count} = [cfg.cond{1} ' > ' cfg.cond{2}];
                    else
                        dirn{count} = [cfg.cond{2} ' > ' cfg.cond{1}];
                    end
                else
                    diff(count) = NaN; allmean1(count) = NaN; allmean2(count) = NaN; clustmean1(count) = NaN; clustmean2(count) = NaN;
                    dirn(count) = {''};
                    warning('No data provided, only summarizing p values')
                end
                count = count+1;
            end
        end
    end
    if isfield(stats{c}.cluster,'negclusters')
        for cc = 1:length(stats{c}.cluster.negclusters)
            if stats{c}.cluster.negclusters(cc).prob < cfg.thresh
                try
                    measure{count} = func2str(data{1}.meas{c});
                catch
                    measure{count} = num2str(c);
                end
                pvalue(count) = stats{c}.cluster.negclusters(cc).prob;
                tmp = stats{c}.cluster.label(find(stats{c}.cluster.negclusterslabelmat == cc));
                tmp = cellcat(', ',tmp,'',1);
                labels{count} = [tmp{:}];
                labels{count}(end-1:end) = [];
                
                if ~isempty(data)
                    allmean1(count) = cfg.test(mean(data{1}.data(:,:,c),2));
                    allmean2(count) = cfg.test(mean(data{2}.data(:,:,c),2));
                    
                    clustmean1(count) = cfg.test(mean(data{1}.data(:,find(stats{c}.cluster.negclusterslabelmat == cc),c),2));
                    clustmean2(count) = cfg.test(mean(data{2}.data(:,find(stats{c}.cluster.negclusterslabelmat == cc),c),2));
                    
                    diff(count) = clustmean1(count)-clustmean2(count);
                    if diff(count) > 0
                        dirn{count} = [cfg.cond{1} ' > ' cfg.cond{2}];
                    else
                        dirn{count} = [cfg.cond{2} ' > ' cfg.cond{1}];
                    end
                else
                    diff(count) = NaN; allmean1(count) = NaN; allmean2(count) = NaN; clustmean1(count) = NaN; clustmean2(count) = NaN;
                    dirn(count) = {''};
                    warning('No data provided, only summarizing p values')
                end
                count = count+1;
            end
        end
    end
    
end

if ~exist('measure','var')
    tbl = [];
else
    
    tbl = table;
    tbl.measure = vert(measure);
    tbl.pvalue = vert(pvalue); %report
    tbl.labels = vert(labels);
    tbl.([cfg.cond{1} '_mean']) = vert(allmean1);
    tbl.([cfg.cond{2} '_mean']) = vert(allmean2);
    tbl.([cfg.cond{1} '_clustermean']) = vert(clustmean1);
    tbl.([cfg.cond{2} '_clustermean']) = vert(clustmean2);
    tbl.cluster_diff = vert(diff);
    tbl.summary = vert(dirn);
end