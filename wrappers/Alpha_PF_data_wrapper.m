function [pfout] = Alpha_PF_data_wrapper(dat,psdrange,searchrange)
% Calculates the alpha peak frequency of the subject using the methods of
% Corcoran et al. (2018). Returns NaNs if no alpha peak detected.
%
% Required inputs:
%   dat is the data in EEGLAB format
% 
% Optional inputs: 
%   psdrange is the frequency range over which to compute the PSD (default 
%       = [1 40])
%   searchrange is the range in which to look for an alpha peak (default =
%       [5 15]) 
%
% NOTE: this function finds one peak frequency for the subject, so it
% returns identical values for all channels. Do NOT use this method with
% a cluster-based permutation test!
%
% Corcoran, A. W., Alday, P. M., Schlesewsky, M., & Bornkessel?Schlesewsky, 
% I. (2018). Toward a reliable, automated method of individual alpha 
% frequency (IAF) quantification. Psychophysiology, 55(7), e13064. 


if ~exist('psdrange','var')
   psdrange = [1 40]; 
end

if ~exist('searchrange','var')
   searchrange = [5 15]; 
end

disp('Computing alpha peak frequency...')
disp('')

[psum] = restingIAF(dat.data,dat.nbchan,3,psdrange,dat.srate,searchrange,11,5);

if ~isnan(psum.paf)
    pf = psum.paf;
else
    pf = psum.cog;
end

pfout = repmat(pf,1,dat.nbchan);
