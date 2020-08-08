function [bp] = Bandpower_data_wrapper(dat,frange,norm_bandpass)
% Computes the power in a given frequency band
%
% Required inputs: 
%   dat: the original data, in EEGLAB format
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

for c = 1:dat.nbchan
   bp(c) = bandpower(dat.data(c,:),dat.srate,frange); 
   if ~strcmpi(norm_bandpass,'no')
      bp(c) = bp(c)/bandpower(dat.data(c,:),dat.srate,norm_bandpass); 
   end
end