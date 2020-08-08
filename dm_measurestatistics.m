function [stats] = ft_measurestatistics(cfg,data)
% ft_measurestatistics does stats on two (eventually or more) outputs of
% ft_applymeasure
%
% Input arguments:
%
% cfg: a config structure with the following options:
%      test: the test you want to use, passed in as a string. Inputs
%         include 'ttest' (for paired samples), 'ttest2' (for unpaired),
%         'signrank', 'ranksum', 'anova', 'rmanova', 'kruskalwallis',
%         'friedman', 'tost', and 'empirical'. Use empirical for comparing one
%         outputs structure with a fourth dimension (either surrogates or
%         subsamples) with another outputs structure with only three
%         dimensions. (default = 'ranksum' or 'kruskal-wallis' depending on
%         the length of data)
%      channel: the channels you want to test - must be an input that works
%         with ft_channelselection. (default = 'all');
%      multcompare: your method of multiple comparison correction.
%         'cluster', 'fdr', 'mean', or 'none' (default = 'cluster')
%      effectsize: do effect size statistics - inputs can be any of the
%         inputs to the 'mes.m' function from the Measures of Effect Size
%         Toolbox (Hentschke and St�ttgen, 2011) (default set based on
%         cfg.test)
%      subs: a vector or cell array of subject indices to choose from. This
%         allows selection of common subjects if some are missing from one 
%         condition (default = all subjects for each element of 'data')
%      cluster: if you chose 'cluster' as the multiple comparison
%         correction method, you can add optional inputs for the
%         permutation test
%              nrand: number of permutations (default = 10000)
%              minnbchan: minimum number of channels in a cluster (default
%              = 1)
%              statfun: the fieldtrip function to use for the permutation
%              test (default = set based on cfg.test)
%              datasetinfo: the channel information used by EasyClusterCorrect
%              (default = from data)
%      eqinterval: equivalence interval for non-inferiority testing. 
% data: a 1 X N cell array, each cell containing one "outputs" structure
%      from ft_applymeasure. Each outputs structure must have the same
%      dimensions
%
%
% Outputs:
%
% stats: a 1 X number of measures cell array, each containing a stats
%    structure with fields depending on the input configuration

%% Set defaults

ft_defaults

if ~cfgcheck(cfg,'test')
    if length(data) == 2
        cfg.test = 'ranksum';
    else
        cfg.test = 'kruskalwallis';
    end
end

if ~cfgcheck(cfg,'channel')
    cfg.channel = {'all'};
end

cfg = setdefault(cfg,'meas',1:length(data{1}.meas));

if ~cfgcheck(cfg,'multcompare')
    cfg.multcompare = 'cluster';
    cfg.cluster = struct;
end

if cfgcheck(cfg,'multcompare','cluster') && ~isfield(cfg,'cluster')
    cfg.cluster = struct;
end

if ~cfgcheck(cfg,'effectsize')
    switch cfg.test
        case {'ttest','ttest2','tost'}
            cfg.effectsize = 'hedgesg';
        case 'ranksum'
            cfg.effectsize = 'auroc';
        case 'signrank'
            cfg.effectsize = 'rbcorr';
        case {'anova','kruskalwallis','friedman'}
            cfg.effectsize = 'psi';
    end
end

switch cfg.test
    case {'ttest','signrank','friedman'}
        mesargsin = {'isDep',1};
    otherwise
        mesargsin = {};
end

if cfgcheck(cfg,'multcompare','cluster') && ~cfgcheck(cfg.cluster,'statfun')
    switch cfg.test
        case 'ttest'
            cfg.cluster.statfun = 'ft_statfun_depsamplesT';
        case 'ttest2'
            cfg.cluster.statfun = 'ft_statfun_indepsamplesT';
        case 'ranksum'
            cfg.cluster.statfun = 'ft_statfun_ranksum';
        case 'signrank'
            cfg.cluster.statfun = 'ft_statfun_fast_signrank';
        case 'anova'
            cfg.cluster.statfun = 'ft_statfun_indepsamplesF'; 
        case 'rmanova'
            cfg.cluster.statfun = 'ft_statfun_depsamplesF';
        case 'kruskalwallis'
            cfg.cluster.statfun = 'ft_statfun_kruskal';
        case 'friedman'
            cfg.cluster.statfun = 'ft_statfun_friedman';
	case 'tost'
	    cfg.cluster.statfun = 'ft_statfun_tost';
    end
end

if cfgcheck(cfg,'multcompare','cluster') && ~cfgcheck(cfg.cluster,'nrand')
    cfg.cluster.nrand = 10000;
end

if cfgcheck(cfg,'multcompare','cluster') && ~cfgcheck(cfg.cluster,'minnbchan')
    cfg.cluster.minnbchan = 1;
end

if isfield(cfg,'eqinterval')
   cfg.cluster.eqinterval = cfg.eqinterval; 
end

cfg = setdefault(cfg,'subs',repmat({'all'},1,length(data)));

%% Select channels
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

if cfgcheck(cfg,'multcompare','cluster')
    cfg.cluster.channel = data{1}.chan(chans);
end

% for c = 1:length(data)
%     % for now, chans is always the second dimension, so no need to figure
%     % out how to work with dimord
%     
%     %dimns = tokenize(data{1}.dimord,'_');
%     %chans = find(strcmpi(dimns,'chan'));
%     data{c}.data = data{c}.data(:,chans,:,:);
% end

%% Select subjects

for c = 1:length(data)
   if ~strcmpi(cfg.subs{c},'all')
       data{c}.data = data{c}.data(cfg.subs{c},:,:,:);
       if isfield(data{c},'sub')
            data{c}.sub = data{c}.sub(cfg.subs{c});
       end
   end
end

%% Calculate stats
for i = cfg.meas
    %dimns = tokenize(data{1}.dimord,'_');
    if ~cfgcheck(cfg,'multcompare','mean')
        for c = 1:length(chans)
            switch cfg.test
                case 'ttest'
                    [~,stats{i}.p(c)] =  ttest(data{1}.data(:,chans(c),i)-data{2}.data(:,chans(c),i));
                case 'ttest2'
                    [~,stats{i}.p(c)] = ttest2(data{1}.data(:,chans(c),i),data{2}.data(:,chans(c),i));
                case 'ranksum'
                    stats{i}.p(c) = ranksum(data{1}.data(:,chans(c),i),data{2}.data(:,chans(c),i));
                case 'signrank'
                    stats{i}.p(c) = signrank(data{1}.data(:,chans(c),i),data{2}.data(:,chans(c),i));
                case 'anova'
                    dat = []; design = [];
                    for cc = 1:length(data)
                        dat = cat(2,dat,horz(data{cc}.data(:,chans(c),i)));
                        design = cat(2,design,ones(1,length(data{cc}.data(:,chans(c),i)))*cc);
                    end
                    stats{i}.p(c) = anovan(dat,design,'off');
                case 'rmanova'
                    % not implemented yet
                case 'kruskalwallis'
                    dat = []; design = [];
                    for cc = 1:length(data)
                        dat = cat(2,dat,horz(data{cc}.data(:,chans(c),i)));
                        design = cat(2,design,ones(1,length(data{cc}.data(:,chans(c),i)))*cc);
                    end
                    stats{i}.p(c) = kruskalwallis(dat,design,'off');
                case 'friedman'
                    dat = [];
                    for cc = 1:length(data)
                        dat(cc,:) = vert(data{cc}.data(:,chans(c),i));
                    end
                    stats{i}.p(c) = friedman(dat,1,'off');
                case 'empirical'
                    % finish later
                case 'tost'
                    [p1,p2] = TOST(data{1}.data(:,chans(c),i),data{2}.data(:,chans(c),i),...
                        cfg.eqinterval(1),cfg.eqinterval(2),0.05);
                    stats{i}.p(c) = max(p1,p2);
            end
            
            switch cfg.test
                case {'ttest','ttest2','ranksum','signrank'}
                    stats{i}.effsize{c} = mes(data{1}.data(:,chans(c),i),data{2}.data(:,chans(c),i),cfg.effectsize,mesargsin{:});
                case {'anova','rmanova','kruskalwallis','friedman'}
                    % not implemented yet
            end
        end
        
        switch cfg.multcompare
            case 'cluster'
                if ~isfield(cfg.cluster,'datasetinfo')
                    if isfield(data{1},'elec')
                        datasetinfo.elec = data{1}.elec;
                    elseif isfield(data{1},'grad')
                        datasetinfo.grad = data{1}.grad;
                    end
                    datasetinfo.label = data{1}.chan;
                else
                    datasetinfo = cfg.cluster.datasetinfo;
                end
                
                for c = 1:length(data)
                    input{c} = data{c}.data(:,:,i)';
                end
                stats{i}.cluster = EasyClusterCorrect(input,datasetinfo,cfg.cluster.statfun,cfg.cluster);
                stats{i}.mask = stats{i}.cluster.mask;
            case 'fdr'
                stats{i}.fdr = fdr(stats{i}.p);
                stats{i}.mask = stats{i}.fdr < 0.05;
        end
    elseif cfgcheck(cfg,'multcompare','mean')
        switch cfg.test
            case 'ttest'
                [~,stats{i}.p] =  ttest(nanmean(data{1}.data(:,:,i),2)-nanmean(data{2}.data(:,:,i),2));
            case 'ttest2'
                [~,stats{i}.p] = ttest2(nanmean(data{1}.data(:,:,i),2),nanmean(data{2}.data(:,:,i),2));
            case 'ranksum'
                stats{i}.p = ranksum(nanmean(data{1}.data(:,:,i),2),nanmean(data{2}.data(:,:,i),2));
            case 'signrank'
                stats{i}.p = signrank(nanmean(data{1}.data(:,:,i),2),nanmean(data{2}.data(:,:,i),2));
            case 'anova'
                dat = []; design = [];
                for cc = 1:length(data)
                    dat = cat(2,dat,horz(nanmean(data{cc}.data(:,:,i),2)));
                    design = cat(2,design,ones(1,length(data{cc}.data(:,1,i)))*cc);
                end
                stats{i}.p = anovan(dat,design,'off');
            case 'rmanova'
                % not implemented yet
            case 'kruskalwallis'
                dat = []; design = [];
                for cc = 1:length(data)
                    dat = cat(2,dat,horz(nanmean(data{cc}.data(:,:,i),2)));
                    design = cat(2,design,ones(1,length(data{cc}.data(:,1,i)))*cc);
                end
                stats{i}.p = kruskalwallis(dat,design,'off');
            case 'friedman'
                dat = [];
                for cc = 1:length(data)
                    dat(cc,:) = vert(nanmean(data{cc}.data(:,:,i),2));
                end
                stats{i}.p = friedman(dat,1,'off');
            case 'empirical'
                % finish later
                
            case 'tost'
                [p1,p2] = TOST(nanmean(data{1}.data(:,:,i),2),nanmean(data{2}.data(:,:,i),2),...
                    cfg.eqinterval(1),cfg.eqinterval(2),0.05);
                stats{i}.p = max(p1,p2);
        end
        switch cfg.test
            case {'ttest','ttest2','ranksum','signrank''tost'}
                stats{i}.effsize = mes(mean(data{1}.data(:,:,i),2),mean(data{2}.data(:,:,i),2),cfg.effectsize,mesargsin{:});
            case {'anova','rmanova','kruskalwallis','friedman'}
                % not implemented yet
        end
    end
    stats{i}.cfg = cfg;
end


