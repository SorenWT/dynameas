function [osciout] = osci_fooof_wrapper(dat,frange,method,varargin)
% Calculates the power-law exponent
%
% Required inputs:
%   dat: the outputs from a fooof transform
%   frange: the frequency range of interest
%
% Recommended inputs:
%   method: whether to compute the oscillatory power based on the residual
%      power spectrum after subtracting the aperiodic fit (method =
%      'resid'), or whether to compute it based on the fit oscillatory
%      peaks (method = 'peaks')
%
% Optional inputs: optional inputs for this function are in the form of
%   name-value pairs


if nargin < 4
    argsin = {[]};
else
    argsin = varargin;
end

if nargin < 3
    method = 'resid';
end

intout = zeros(1,length(dat));

disp(' ')
ft_progress('init','text','Computing oscillatory power...')


f = dat(1).f;
findx = intersect(find(f > frange(1)),find(f < frange(2)));

for i = 1:length(dat)
    switch method
        case 'resid'
            osciout(i) = simps(f(findx),dat(i).resid_pxx(findx))/numel(findx);
        case 'peaks'
    end
end

ft_progress('close')