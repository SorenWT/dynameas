function [degreeout] = degree_connectivity_wrapper(dat,frange,varargin)
% Calculates the power-law exponent, as in Northoff et al. (2020)
% 
% Required inputs:
%   dat: the outputs from a connectivity transform
% 
% Recommended inputs:
%   frange: the frequency range over which to average connectivity (default =
%       whole range of the data)
% 
% Optional inputs: optional inputs for this function are in the form of
%   name-value pairs
%   'pmethod', 'periodogram': uses periodogram-based PSD estimation instead 
%       of Welch's method
%   'interp','off': turns off frequency interpolation


if nargin < 3
    argsin = {[]};
else
    argsin = varargin;
end

if nargin < 2
    frange = [min(dat.freq)-1 max(dat.freq)+1];
end

dat.nbchan = length(dat.label);

degreeout = zeros(1,dat.nbchan);

disp(' ')
ft_progress('init','text','Computing node degree...')

if ischar(frange)
    frange = eval(frange);
end

findx = find(dat.freq>=frange(1) & dat.freq<=frange(2));


connparam = [dat.cfg.method 'spctrm'];

%for c = 1:dat.nbchan
        %ft_progress(c/dat.nbchan,'Processing channel %d out of %d',c,dat.nbchan);
    degreeout = sum(mean(dat.(connparam)(:,:,findx),3),1);
        %PLEout(c) = JF_power_law(dat.data(c,:),dat.srate,frange(1),frange(2),argsin{:});
%end
ft_progress('close')