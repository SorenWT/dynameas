function freqdata = transform_timefreq(cfg,data)

cfg = setdefault(cfg,'method','mtmconvol');
cfg = setdefault(cfg,'load','yes');
cfg = setdefault(cfg,'save','no');
cfg = setdefault(cfg,'fldr',pwd);
cfg = setdefault(cfg,'pad','nextpow2');
cfg = setdefault(cfg,'keeptrials','yes'); 
cfg = setdefault(cfg,'output','fourier');
cfg = setdefault(cfg,'computepaf','no');

disp(' ')
disp('Doing time-frequency transformation...')
fname = data.filename;

ftdata = eeglab2fieldtrip(data,'preprocessing','none');
ftdata = ft_concat(ftdata);


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
        %if min(cfg.foi) >= 30
           cfg = setdefault(cfg,'taper','dpss'); 
           cfg = setdefault(cfg,'tapsmofrq',0.2*cfg.foi);
           cfg = setdefault(cfg,'t_ftimwin',5./cfg.foi);
        %else
        %    cfg = setdefault(cfg,'taper','hanning');
        %    cfg = setdefault(cfg,'t_ftimwin',repmat(3/min(cfg.foi),1,length(cfg.foi)));
        %end
        % 50% overlap on the lowest frequency
        cfg = setdefault(cfg,'toi',linspace(data.times(1),data.times(end),2*data.times(end)./max(cfg.t_ftimwin)));
    case 'wavelet'
        cfg = setdefault(cfg,'foi',2:2:50);
        cfg = setdefault(cfg,'toi',linspace(data.times(1),data.times(end),round(length(data.times)/25)));
        cfg = setdefault(cfg,'width',3);
    case 'irasa'
        % do this later
end

if exist(fullfile(cfg.fldr,[fname '_timefreq_' cfg.method '.mat']),'file') && strcmpi(cfg.load,'yes')
    freqdata = parload(fullfile(cfg.fldr,[fname '_timefreq_' cfg.method '.mat']),'freqdata');
    cfg.save = 'no'; % no need to resave if we just loaded
else
    tmpcfg = cfg; tmpcfg = rmfield(cfg,{'load','save'}); 
    freqdata = ft_freqanalysis(tmpcfg,ftdata);
end

freqdata.cumtapcnt = ones(size(freqdata.fourierspctrm,4),1);

if strcmpi(cfg.computepaf,'yes')
    [~,psum] = Alpha_PF_data_wrapper(data);
    
    freqdata.cfg.dm.pfinfo = psum;
    freqdata.cfg.dm.alphacog = psum.cog;
    freqdata.cfg.dm.alphapf = psum.paf;
    if isfield(psum,'iaw') && ~any(isnan(psum.iaw))
    freqdata.cfg.dm.iaw = psum.f(psum.iaw);
    else
        freqdata.cfg.dm.iaw = [8 15];
    end
end


if cfgcheck(cfg,'save','yes') 
    save(fullfile(cfg.fldr,[data.filename '_timefreq_' cfg.method '.mat']),'freqdata','-v7.3')
end


end