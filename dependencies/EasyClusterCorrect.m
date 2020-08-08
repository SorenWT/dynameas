function [stats] = EasyClusterCorrect(data,datasetinfo,statfun,opts)
% EasyClusterCorrect implements fieldtrip's cluster correction algorithm
% across channels only, for an arbitrary measurement
%
% Input arguments:
%      data: a 1xn cell array containing the data for testing. Each cell
%         should be a matrix of the form channels x observations, or a 3-D
%         matrix of the form channels x observations x time (for
%         compatibility with the 2-D form)
%      datasetinfo: a structure containing the relevant information for
%         clustering. For the following data types, this should be:
%         eeg: requires fields 'elec' (a fieldtrip electrode structure) and
%            'label'
%         meg: requires fields 'grad' (a fieldtrip gradiometer structure)
%            and 'label'
%         source: currently only region-based correction supported.
%            Requires fields 'atlas' (a fieldtrip atlas structure) and
%            'atlasname' (currently only 'aal','mmp', and 'yeo' suppported
%            - for 'yeo' the atlas must be modified to include 'pos' and
%            'tri' fields from an HCP subject's source model)
%      statfun: the fieldtrip statfun used for the statistics (e.g.
%         ft_statfun_signrank)
%      opts: optional argument which changes some parameters of the
%         clustering. Can have the following fields:
%         nrand: the number of permutations to carry out (default = 2000)
%         minnbchan: the minimum number of neighbours to add a channel to a
%            cluster (default = 1)
%
% Outputs:
%      stats: a fieldtrip stats structure returned from
%         ft_timelockstatistics


if nargin < 4
    opts = struct;
end

opts = setdefault(opts,'nrand',2000);

opts = setdefault(opts,'minnbchan',1);

opts = setdefault(opts,'type','Spearman');

opts = setdefault(opts,'channel',{'all'});

if length(datasetinfo.label) >= 32
    opts = setdefault(opts,'distmethod','distance');
elseif length(datasetinfo.label) > 1
    opts = setdefault(opts,'distmethod','triangulation');
end

if ~isfield(opts,'designtype')
   switch statfun
       case {'ft_statfun_signrank','ft_statfun_fast_signrank','ft_statfun_spearman',...
               'ft_statfun_depsamplesT','ft_statfun_friedman','ft_statfun_partspearman',...
               'ft_statfun_depsamplesregrT','ft_statfun_actvsblT'}
            opts.designtype = 'dep';
       otherwise
           opts.designtype = 'indep';
   end
               
end

fields = fieldnames(datasetinfo);
badlabels = [];

for c = 1:length(data)
    tlock{c} = struct;
    if size(data{c},3) == 1
        tlock{c}.avg = [mean(data{c},2) mean(data{c},2) mean(data{c},2)];
        tlock{c}.trial = cat(3,data{c}',data{c}',data{c}');
        tlock{c}.time = [-1 0 1];
    else
        tlock{c}.avg = mean(data{c},2);
        tlock{c}.trial = t3d(data{c});
        tlock{c}.time = linspace(0,0.5,size(data{c},3));
    end
    tlock{c}.dimord = 'rpt_chan_time';
    
    
    if isfield(datasetinfo,'elec') || isfield(datasetinfo,'grad')
        for q = 1:length(fields)
            tlock{c}.(fields{q}) = datasetinfo.(fields{q});
        end
    elseif isfield(datasetinfo,'atlasname')
        switch datasetinfo.atlasname
            case 'aal'
                tlock{c}.label = datasetinfo.atlas.tissuelabel;
            case {'mmp','yeo'}
                tlock{c}.label = datasetinfo.atlas.parcellationlabel;
        end
    end
    for cc = 1:length(tlock{c}.label)
        if isempty(find(~isnan(tlock{c}.trial(:,cc,1))))
            badlabels = [badlabels {tlock{c}.label{cc}}];
        end
    end
    %tlock{c}.label = datasetinfo.label;
    %tlock{c}.grad = datasetinfo.grad;
    tlock{c}.fsample = 1;
    tlock{c}.sampleinfo = [[1:size(tlock{c}.trial,1)]' [1:size(tlock{c}.trial,1)]'+3];
end

if isfield(datasetinfo,'grad') || isfield(datasetinfo,'elec')
    cfg = []; cfg.method = opts.distmethod;
    if length(datasetinfo.label) > 1
        neighbs = ft_prepare_neighbours(cfg,datasetinfo);
    end
else
    switch datasetinfo.atlasname
        case 'aal'
            tissue = datasetinfo.atlas.tissue;
            tissue(find(tissue == 0)) = NaN;
            
            neighbs = struct;
            for c = 1:length(datasetinfo.atlas.tissuelabel)
                neighbs(c).label = datasetinfo.atlas.tissuelabel{c};
                neighbs(c).neighblabel = cell(1,1);
            end
            
            for c = 1:3
                edges{c} = diff(tissue,1,c);
                edges{c} = edges{c} > 0;
                boundaryindx = find(edges{c});
                for cc = 1:length(boundaryindx)
                    [i1,i2,i3] = ind2sub(size(edges{c}),boundaryindx(cc));
                    indx = {i1 i2 i3};
                    indx{c} = indx{c}+1;
                    roi1val = tissue(i1,i2,i3);
                    roi2val = tissue(indx{:});
                    neighbs(roi1val).neighblabel = [neighbs(roi1val).neighblabel neighbs(roi2val).label];
                    neighbs(roi2val).neighblabel = [neighbs(roi2val).neighblabel neighbs(roi1val).label];
                end
            end
            
            for c = 1:length(neighbs)
                neighbs(c).neighblabel(1) = [];
                neighbs(c).neighblabel = unique(neighbs(c).neighblabel);
            end
        case {'mmp','yeo'}
            
            if strcmpi(datasetinfo.atlasname,'yeo')
                datasetinfo.atlas.parcellationlabel = cellstr(num2str([1:8004]'));
                datasetinfo.atlas.parcellation = datasetinfo.atlas.parcels;
            end
            
            pos = datasetinfo.atlas.pos;
            tri = datasetinfo.atlas.tri;
            
            vox_neighbs = cell(1,length(pos));
            
            for c = 1:length(pos)
                inds = find(tri == c);
                for cc = 1:length(inds)
                    [i1,i2] = ind2sub(size(tri),inds(cc));
                    vox_neighbs{c} = [vox_neighbs{c} datasetinfo.atlas.tri(i1,except(1:3,i2))];
                end
                vox_neighbs{c} = unique(vox_neighbs{c});
            end
            
            reg_neighbs = cell(1,length(datasetinfo.atlas.parcellationlabel));
            for c = 1:length(datasetinfo.atlas.parcellationlabel)
                reg_neighbs{c} = cat(2,vox_neighbs{find(datasetinfo.atlas.parcellation == c)});
                for cc = 1:length(reg_neighbs{c})
                    reg_neighbs{c}(cc) = datasetinfo.atlas.parcellation(reg_neighbs{c}(cc));
                end
                reg_neighbs{c} = unique(reg_neighbs{c});
                reg_neighbs{c}(find(reg_neighbs{c} == c)) = [];
            end
            
            neighbs = struct;
            for c = 1:length(reg_neighbs)
                neighbs(c).label = datasetinfo.atlas.parcellationlabel{c};
                neighbs(c).neighblabel = {datasetinfo.atlas.parcellationlabel{reg_neighbs{c}}};
            end
            
            vox_neighbs = [];
            reg_neighbs = [];
            pos = [];
            tri = [];
    end
end

badlabels = unique(badlabels);
for c = 1:length(data)
    cfg = []; cfg.channel = except(1:length(tlock{c}.label),Subset_index(tlock{c}.label,badlabels));
    tlock{c} = ft_selectdata(cfg,tlock{c});
end

if length(datasetinfo.label) > 1
    rmneighbs = zeros(1,length(neighbs));
    for c = 1:length(neighbs)
        if ~isempty(find(strcmpi(neighbs(c).label,badlabels)))
            rmneighbs(c) = 1;
        end
    end
    neighbs(find(rmneighbs)) = [];
end

cfg = []; cfg.method = 'montecarlo'; cfg.statistic = statfun;
cfg.correctm = 'cluster'; cfg.clusterstatistic = 'maxsum';
if length(datasetinfo.label) > 1
    cfg.neighbours = neighbs;
end

if strcmpi(statfun,'ft_statfun_indepsamplesF')
    cfg.tail = 1; cfg.clustertail = 1; cfg.alpha = 0.025; cfg.clusteralpha = 0.05;
else
    cfg.tail = 0; cfg.clustertail = 0; cfg.alpha = 0.025; cfg.clusteralpha = 0.05;
end
cfg.numrandomization = opts.nrand; cfg.spmversion = 'spm12'; cfg.minnbchan = opts.minnbchan;

cfg.type = opts.type;

if isfield(opts,'parpool')
    cfg.parpool = opts.parpool;
end

if ~isfield(opts,'external')
    design = zeros(1,size(cat(2,data{:}),2));
    currindx = 1;
    for c = 1:length(data)
        design(1,currindx:size(data{c},2)) = 1;
        currindx = currindx+size(data{c},2);
    end
    design = design+1;
else
    design = opts.external;
end
%design(1,(size(data{1},2)+1):(size(data{1},2) + size(data{2},2)))= 2;
if strcmpi(opts.designtype,'indep')
    cfg.design = design; cfg.ivar = 1;
else
    siz = cellfun(@(dat)size(dat,2),data,'UniformOutput',true);
    des = [];
    for i = 1:length(siz)
       des = cat(2,des,1:siz(i)); 
    end
    cfg.design(1,:) = des;
    cfg.design(2,:) = design; 
    cfg.uvar = 1; cfg.ivar = 2;
end

for c = 1:length(tlock)
    newdata{c} = tlock{c}.trial(:,:,1)';
end

if isfield(opts,'partial')
   cfg.partial = opts.partial; 
end

%if any(cell2mat(cellfun(@any,cellfun(@isnan,newdata,'UniformOutput',false),'UniformOutput',false)))
%    cfg.nanrand = 'yes';
%end

if isfield(opts,'eqinterval')
   cfg.eqinterval = opts.eqinterval; 
end

%if isfield(datasetinfo,'grad')
%    cfg.channel = {'MEG'};
%else
%    cfg.channel = {'EEG'};
%end
cfg.channel = opts.channel;

cfg.latency = [0 0.5];
cfg.parameter = 'trial';

stats = ft_timelockstatistics(cfg,tlock{:});

% fixing outputs to be more intuitive

if isfield(stats,'posclusters')
    for i = 1:length(stats.posclusters)
        stats.posclusters(i).prob = 2*stats.posclusters(i).prob;
        if stats.posclusters(i).prob > 1
           stats.posclusters(i).prob = 1; 
        end
    end
end

if isfield(stats,'negclusters')
    for i = 1:length(stats.negclusters)
        stats.negclusters(i).prob = 2*stats.negclusters(i).prob;
        if stats.negclusters(i).prob > 1
           stats.negclusters(i).prob = 1; 
        end
    end
end

stats.prob = 2*stats.prob; stats.prob(stats.prob > 1) = 1;
stats.message = 'p values multiplied by 2 to reflect overall two-tailed alpha level of 0.05';