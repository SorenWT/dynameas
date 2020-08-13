function [PLEout] = lowpsdwe_data_wrapper(dat,frange,varargin)
% Calculates the power-law exponent using the enhanced periodogram method
% of Eke et al. (2000)
% 
% Required inputs:
%   dat: the raw data in EEGLAB format
% 
% Recommended inputs:
%   frange: the frequency range over which to calculate the PLE (default =
%       [1 50])
% 
% Optional inputs: optional inputs for this function are in the form of
%   name-value pairs
%   'pmethod', 'periodogram': uses periodogram-based PSD estimation instead 
%       of Welch's method
%   'interp','off': turns off frequency interpolation



if ~exist('frange','var')
   frange = [1 50];
end

if nargin < 3
    argsin = {[]};
else
    argsin = varargin;
end


PLEout = zeros(1,dat.nbchan);

disp(' ')
ft_progress('init','text','Computing power law exponent...')

for c = 1:dat.nbchan
        ft_progress(c/dat.nbchan,'Processing channel %d out of %d',c,dat.nbchan);
    fprintf([num2str(c) ' ']);
    PLEout(c) = lowpsdwe(dat.data(c,:),dat.srate,frange(1),frange(2),argsin{:});
end
ft_progress('close')