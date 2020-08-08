function [bpout] = Alphapower_individ_data_wrapper(dat,norm_bandpass,psdrange,searchrange)
% This function computes individualized alpha power based on the peak
% frequency and width determined from the method described in Corcoran et
% al. (2018). 
%
% Required inputs: 
%   dat: the raw data, in EEGLAB format
%
% Recommended inputs: 
%   norm_bandpass: determines whether to normalize the power by the total
%       power within some bandwidth. Input a tuple to specify this range,
%       input 'no' to use absolute power (default = 'no')
%
% Optional inputs: 
%   psdrange: the frequency range over which to compute the PSD (default =
%       [1 40])
%   searchrange: the range in which to look for an alpha peak (default =
%       [5 15]) 


% Corcoran, A. W., Alday, P. M., Schlesewsky, M., & Bornkessel?Schlesewsky, 
% I. (2018). Toward a reliable, automated method of individual alpha 
% frequency (IAF) quantification. Psychophysiology, 55(7), e13064.



if ~exist('norm_bandpass','var')
   norm_bandpass = 'no'; 
end

if ~exist('psdrange','var')
   psdrange = [1 40]; 
end

if ~exist('searchrange','var')
   searchrange = [5 15]; 
end

disp('Computing individualized alpha power...')
disp('')

[psum,~,f] = restingIAF(dat.data,dat.nbchan,3,psdrange,dat.srate,searchrange,11,5);

if ~isnan(psum.paf)
    pf = psum.paf;
else
    pf = psum.cog;
end

try
    iarange = [f(psum.iaw(1)) f(psum.iaw(2))];
    bpout = Bandpower_EEG_wrapper(dat,iarange,norm_bandpass);
catch
    bpout = NaN(1,dat.nbchan);
end

