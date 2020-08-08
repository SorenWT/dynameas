function teresults = TE_estimation(data,lags,highpass)
% Inputs:
%     data: a fieldtrip data structure, epoched
%     lags: a scalar or vector of lags (in samples) at which to compute
%        transfer entropy (default = 1:ceil(fs))
%     highpass: the highpass frequency of the data (used for
%     cfg.actthrvalue)
%
% Outputs:
%     teresults: transfer entropy results. Important fields are:
%         TEmat: the actual transfer entropy values. The format is channels
%         x music pieces, or channels x shifts for single-trial data
%         MImat: mutual information values. Format is the same as TEmat
%         (channels x music pieces/shifts)


if ~exist('lags','var')
    lags = 1:ceil(fs);
end

if ~exist('highpass','var')
   highpass = 1; 
end


for i = 1:length(lags)
    cfg = [];
    cfg.ensemblemethod = 'no';
    cfg.channel = data.label;
    cfg.toi = [-Inf Inf];
    cfg.ragtaurange = [0.2 0.5]; cfg.ragdim = 2:12;
    cfg.actthrvalue = (1/highpass)*data.fsample;
    cfg.minnrtrials = 1; cfg.maxlag = ceil(0.5*size(data.trial{1},2)); cfg.repPred = ceil(0.2*size(data.trial{1},2));
    cfg.predicttime_u = 1000*lags(i)/data.fsample;
    cfg.flagNei = 'Mass'; cfg.sizeNei = 4;
    
    % remove electrodes with outlier ACT
%     act = TEgetACT(cfg,data);
%     act_thr = max(mean(mean(act,3),1)+3*std(mean(act,3),[],1));
%     for ii = 1:size(act,1)
%         goodact(ii) = mean(act(ii,2,:),3) < act_thr;
%     end
%     tmpcfg = []; tmpcfg.channel = [vert(labels(goodact)); {'Music'}];
%     seldata = ft_selectdata(tmpcfg,data);
%     cfg.actthrvalue = act_thr;
%     cfg.sgncmb = cat(2,repmat({'Music'},length(seldata.label)-1,1),vert(seldata.label(1:end-1)));
    
    %dataprep = TEprepare(cfg,seldata);
    dataprep = TEprepare(cfg,data);
    
    cfg = []; cfg.optdimusage = 'maxdim';
    cfg.fileidout = 'test';
    cfg.surrogatetype = 'trialshuffling';
    cfg.shifttest = 'yes';
    cfg.numpermutation = 1000;
    
    teresults{i} = TEsurrogatestats(cfg,dataprep);
    
%     teresults{i}.act = act;
%     teresults{i}.sgncmb = cat(2,repmat({'Music'},length(labels),1),vert(labels));
%     %teresults{i}.trials(125,:) = teresults.trials(124,:);
%     tmp = NaN(length(goodact),size(teresults{i}.TEmat,2));
%     tmp(goodact,:) = teresults{i}.TEmat;
%     teresults{i}.TEmat = tmp;
%     tmp = NaN(length(goodact),size(teresults{i}.TEmat,2));
%     tmp(goodact,:) = teresults{i}.MImat;
%     teresults{i}.MImat = tmp;
end





