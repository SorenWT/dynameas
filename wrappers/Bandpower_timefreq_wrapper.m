function [bp] = Bandpower_timefreq_wrapper(dat,frange,norm_bandpass)
% Computes the power in a given frequency band
%
% Required inputs: 
%   dat: the timefreq-transformed data (i.e. Fieldtrip freq format)
%   frange: the frequency range in which you want to compute power, as a
%       tuple
% 
% Recommended inputs:
%   norm_bandpass: determines whether to normalize the power by the total
%       power within some bandwidth. Input a tuple to specify this range,
%       input 'no' to use absolute power (default = 'no')

if ~exist('norm_bandpass','var')
   norm_bandpass = 'no'; 
end

if ischar(frange)
   frange = eval(frange); 
end

disp(' ')

ft_progress('init','text','Computing power...')

psd = squeeze(nanmean(abs(dat.fourierspctrm),4));

findx = intersect(find(dat.freq > frange(1)),find(dat.freq < frange(2)));

%bp = trapz(dat.freq(findx),psd(:,findx)');
bp = norm(psd(:,findx)'.^2)./numel(findx);

if ~strcmpi(norm_bandpass,'no')
    findx_norm = intersect(find(dat.freq > norm_bandpass(1)),find(dat.freq < norm_bandpass(2)));
   %bp = bp./trapz(dat.freq(findx_norm),psd(:,findx_norm)');
    bp = bp./norm(psd(:,findx_norm)'.^2)./numel(findx_norm);
end

ft_progress('close');