function tbl = ft_measurestatsummary(cfg,data,stats)
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
% data: a cell array of outputs from ft_applymeasure - these should be in
%      the SAME ORDER as input to ft_measurestatistics
% stats: the output from ft_measurestatistics
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

cfg = setdefault(cfg,'thresh',0.1);
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
            if stats{c}.cluster.posclusters(cc).prob < cfg.thresh/2
                measure{count} = func2str(data{1}.meas{c});
                pvalue(count) = stats{c}.cluster.posclusters(cc).prob;
                tmp = stats{c}.cluster.label(find(stats{c}.cluster.posclusterslabelmat == cc));
                tmp = cellcat(', ',tmp,'',1);
                labels{count} = [tmp{:}];
                labels{count}(end-1:end) = [];
                
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
                count = count+1;
            end
        end
    end
    if isfield(stats{c}.cluster,'negclusters')
        for cc = 1:length(stats{c}.cluster.negclusters)
            if stats{c}.cluster.negclusters(cc).prob < cfg.thresh/2
                measure{count} = func2str(data{1}.meas{c});
                pvalue(count) = stats{c}.cluster.negclusters(cc).prob;
                tmp = stats{c}.cluster.label(find(stats{c}.cluster.negclusterslabelmat == cc));
                tmp = cellcat(', ',tmp,'',1);
                labels{count} = [tmp{:}];
                labels{count}(end-1:end) = [];
                
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
                count = count+1;
            end
        end
    end
    
end


tbl = table;
tbl.measure = vert(measure);
tbl.pvalue = vert(pvalue*2); %report
tbl.labels = vert(labels);
tbl.([cfg.cond{1} '_mean']) = vert(allmean1);
tbl.([cfg.cond{2} '_mean']) = vert(allmean2);
tbl.([cfg.cond{1} '_clustermean']) = vert(clustmean1);
tbl.([cfg.cond{2} '_clustermean']) = vert(clustmean2);
tbl.cluster_diff = vert(diff);
tbl.summary = vert(dirn);