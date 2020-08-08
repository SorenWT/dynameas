function [filtdata] = transform_filter(cfg,data)    
% applies band-pass filtering to the data
% 
% cfg options
% 
% Required fields: 
%   bandpass: a tuple specifying the frequency range to bandpass the data
%       into
%
% Recommended fields: none
% 
% Optional fields: 
%   save: save the filtered data to file (default = 'yes')
%   load: if present, load the filtered data from file (default = 'yes')


cfg = setdefault(cfg,'save','no');
cfg = setdefault(cfg,'load','yes');

if exist(fullfile(cfg.fldr,[fname '_filt.mat']),'file') && strcmpi(cfg.load,'yes')
    filtdata = parload(fullfile(cfg.fldr,[fname '_filt.mat']),'irasadata');
else
    data = eeglab2fieldtrip(data,'preprocessing','none');
    
    tmpcfg = []; tmpcfg.bpfilter = 'yes'; tmpcfg.bpfreq = cfg.bandpass;
    tmpcfg.bpinstabilityfix = 'split';
    
    filtdata = ft_preprocessing(tmpcfg,data);
    
    filtdata = ft2eeglab(filtdata); 
end

if strcmpi(cfg.save,'yes')
    save(fullfile(cfg.fldr,[fname '_filt.mat']),'filtdata','-v7.3');
end
end