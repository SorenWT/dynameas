function [seedconnout] = seedbased_connectivity_wrapper(dat,frange,seedchans,varargin)
% Calculates seed-based connectivity from a given channel or set of
% channels to others
% 
% Required inputs:
%   dat: the outputs from a connectivity transform
%   seedchans: the channel or channels to use as a seed, specified in any
%       way that fieldtrip understands (e.g. labels, indices, etc)
% 
% Recommended inputs:
%   frange: the frequency range over which to average connectivity (default =
%       whole range of the data)
% 


if nargin < 3
    argsin = {[]};
else
    argsin = varargin;
end

if nargin < 2
    frange = [min(dat.freq)-1 max(dat.freq)+1];
end

dat.nbchan = length(dat.label);

seedconnout = zeros(1,dat.nbchan);

disp(' ')
ft_progress('init','text','Computing seed-based connectivity...')

if ischar(frange)
   frange = eval(frange); 
end

findx = find(dat.freq>=frange(1) & dat.freq<=frange(2));

connparam = [dat.cfg.method 'spctrm'];

if iscell(seedchans) || ischar(seedchans)
   seedchans = ft_channelselection(seedchans,dat.elec);
   [~,seedchans] = match_str(seedchans,dat.elec.label);
end

seedconnout = mean(mean(dat.(connparam)(seedchans,:,findx),3),1);

ft_progress('close')