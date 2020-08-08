function envdata = transform_envelope(cfg,data)
% takes the amplitude envelope of the data. Generally this is done as a
% second transform following filtering. 
% 
% cfg options:
% 
% Required fields: none
%
% Recommended fields: none
%
% Optional fields: 
%   method: 'fieldtrip' or 'irasa'. IRASA is only used when you
%       specifically want an oscillatory- or fractal-component-specific 
%       amplitude envelope. Otherwise 'fieldtrip' is used. If IRASA is 
%       used, the fields 'winsize(default =
%       'fieldtrip')
%   save: save the enveloped data to file (default = 'yes')
%   load: if present, load the enveloped data from file (default = 'yes')


cfg = setdefault(cfg,'method','fieldtrip');
cfg = setdefault(cfg,'load','yes');
cfg = setdefault(cfg,'save','no');
if cfgcheck(cfg,'method','irasa')
    cfg = setdefault(cfg,'winsize',3);
    cfg = setdefault(cfg,'overlap',0.95);
    cfg = setdefault(cfg,'hset',[1.1:0.05:1.95 2.05:0.05:2.9]);
end

disp(' ')
disp('Getting amplitude envelope...')
fname = data.filename;

if exist(fullfile(cfg.fldr,[fname '_envelope_' cfg.method '.mat']),'file') && strcmpi(cfg.load,'yes')
    envdata = parload(fullfile(cfg.fldr,[fname '_envelope_' cfg.method '.mat']),'envdata');
else
    switch cfg.method
        case 'fieldtrip'
            data = eeglab2fieldtrip(data,'preprocessing','none');
            tmpcfg = []; tmpcfg.hilbert = 'abs';
            envdata = ft_preprocessing(tmpcfg,data);
            envdata = ft2eeglab(envdata);
        case 'irasa'
            [~,envdata] = IRASA_window(data.data,data.srate,'winsize',cfg.winsize,'overlap',cfg.overlap,'hset',cfg.hset);
            envdata(1).srate = 1/(cfg.winsize*(1-cfg.overlap));
    end
end

if cfgcheck(cfg,'save','yes')
    save(fullfile(cfg.fldr,[data.filename '_envelope_' cfg.method '.mat']),'envdata','-v7.3')
end

end