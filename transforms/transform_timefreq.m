function freqdata = transform_timefreq(cfg,data)

cfg = setdefault(cfg,'method','mtmconvol');
cfg = setdefault(cfg,'load','yes');
cfg = setdefault(cfg,'save','no');

switch cfg.method
%     case 'mtmfft'
%         cfg = setdefault(cfg,'foi',1:2:50);
%         if min(cfg.foi) >= 30
%            cfg = setdefault(cfg,'taper','dpss'); 
%            cfg = setdefault(cfg,'tapsmofrq',2);
%         else
%             cfg = setdefault(cfg,'taper','hanning');
%         end
    case 'mtmconvol'
        cfg = setdefault(cfg,'foi',2:2:50);
        if min(cfg.foi) >= 30
           cfg = setdefault(cfg,'taper','dpss'); 
           cfg = setdefault(cfg,'tapsmofrq',2);
        else
            cfg = setdefault(cfg,'taper','hanning');
        end
        cfg = setdefault(cfg,'toi',linspace(data.times(1),data.times(end),round(length(data.times)/25)));
        cfg = setdefault(cfg,'t_ftimwin',repmat(3/min(cfg.foi),1,length(cfg.foi)));
    case 'wavelet'
        cfg = setdefault(cfg,'foi',2:2:50);
        cfg = setdefault(cfg,'toi',linspace(data.times(1),data.times(end),round(length(data.times)/25)));
        cfg = setdefault(cfg,'width',3);
    case 'irasa'
        % do this later
end
cfg = setdefault(cfg,'keeptrials','yes'); 
cfg = setdefault(cfg,'output','fourier');

disp(' ')
disp('Doing time-frequency transformation...')
fname = data.filename;
ftdata = eeglab2fieldtrip(data,'preprocessing','none');

if exist(fullfile(cfg.fldr,[fname '_timefreq_' cfg.method '.mat']),'file') && strcmpi(cfg.load,'yes')
    freqdata = parload(fullfile(cfg.fldr,[fname '_timefreq_' cfg.method '.mat']),'envdata');
else
    tmpcfg = cfg; tmpcfg = rmfield(cfg,{'load','save'}); 
    freqdata = ft_freqanalysis(tmpcfg,ftdata);
end

if cfgcheck(cfg,'save','yes')
    save(fullfile(cfg.fldr,[data.filename '_timefreq_' cfg.method '.mat']),'freqdata','-v7.3')
end

end