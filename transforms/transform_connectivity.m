function conndata = transform_connectivity(cfg,data)

cfg = setdefault(cfg,'method','funccon');
cfg = setdefault(cfg,'load','yes');
cfg = setdefault(cfg,'save','no');
if strcmpi(cfg,'method','transen')
   cfg = setdefault(cfg,'lags',1:round(0.1*data.srate); 
   cfg = setdefault(cfg,'highpass',1);
   cfg = setdefault(cfg,'epochlen',4); 
end

disp(' ')
disp('Creating connectivity matrix...')
fname = data.filename;

if exist(fullfile(cfg.fldr,[fname '_conn_' cfg.method '.mat']),'file') && strcmpi(cfg.load,'yes')
    conndata = parload(fullfile(cfg.fldr,[fname '_conn_' cfg.method '.mat']),'envdata');
else
    ftdata = eeglab2fieldtrip(data,'preprocessing','none');
    switch cfg.method
        case 'transen'
            data = ft_cont2epoch(data,cfg.epochlen);
            warning('Parameters are optimized for M/EEG data. Use with fMRI at your own risk')
            teresults = TE_estimation(data,cfg.lags,cfg.highpass);
            conndata = []; conndata.corr = teresults.TEmat; 
            conndata.label = data.label; conndata.elec = data.elec; conndata.dimord = 'chan_chan';
        case 'funccon'
            conndata = []; conndata.corr = corr(data.trial{1}'); 
            conndata.label = data.label; conndata.elec = data.elec; conndata.dimord = 'chan_chan';
        otherwise
            tmpcfg = cfg; tmpcfg = rmfield(tmpcfg,{'load','save'});
            conndata = ft_connectivityanalysis(tmpcfg,data);
    end
end

if cfgcheck(cfg,'save','yes')
    save(fullfile(cfg.fldr,[data.filename '_conn_' cfg.method '.mat']),'conndata','-v7.3')
end

end