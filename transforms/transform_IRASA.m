function [irasadata] = transform_IRASA(cfg,dat)
% transform_IRASA converts raw data to an IRASA power spectrum structure
% by default, the IRASA decomposition is done in a sliding window
% 
% cfg options: 
% 
% Required fields: none
%
% Recommended fields: none
%
% Optional fields: 
%   winsize: the size (in seconds) of the sliding window for IRASA 
%       decomposition (default = 10 seconds)
%   overlap: the amount of overlap between windows, as a value between 0
%       and 1 (default = 0)
%   hset: the resampling factors used by the IRASA algorithm (default =
%       [1.1:0.05:1.95 2.05:0.05 2.9])
%   save: save the IRASA spectra to file - recommended, as the method takes
%       a long time (default = 'yes')
%   load: if the IRASA spectra have already been computed, load them from
%       file rather than calculating it again (default = 'yes')
%
% Wen, H., & Liu, Z. (2016). Separating Fractal and Oscillatory Components 
% in the Power Spectrum of Neurophysiological Signal. Brain Topography, 
% 29(1), 13?26.



cfg = setdefault(cfg,'winsize',10);
cfg = setdefault(cfg,'overlap',0);
cfg = setdefault(cfg,'hset',[1.1:0.05:1.95 2.05:0.05:2.9]);
cfg = setdefault(cfg,'save','yes');
cfg = setdefault(cfg,'load','yes');

disp(' ')
disp('Performing IRASA...')
fname = dat.filename;

if exist(fullfile(cfg.fldr,[fname '_IRASA_specs.mat']),'file') && strcmpi(cfg.load,'yes')
    irasadata = parload(fullfile(cfg.fldr,[fname '_IRASA_specs.mat']),'irasadata');
else
    irasadata = IRASA_window(dat.data,dat.srate,'winsize',cfg.winsize,'overlap',cfg.overlap,'hset',cfg.hset);
end

if strcmpi(cfg.save,'yes')
    save(fullfile(cfg.fldr,[fname '_IRASA_specs.mat']),'irasadata','-v7.3');
end