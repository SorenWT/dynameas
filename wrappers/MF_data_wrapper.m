function [MFOut] = MF_data_wrapper(dat,frange)
% Computes the median frequency of the data
%
% Required inputs: 
%   dat: the raw data in EEGLAB format
%   
% Recommended inputs:
%   frange: the frequency range over which to compute the median frequency
%       (default = [1 50])

if ~exist('frange','var')
    frange = [1 50];
end


MFOut = zeros(1,dat.nbchan);

disp(' ')
ft_progress('init','text','Computing median frequency...')

for c = 1:dat.nbchan
        ft_progress(c/dat.nbchan,'Processing channel %d out of %d',c,dat.nbchan);
    MFOut(c) = medfreq(dat.data(c,:),dat.srate,frange);
end
ft_progress('close')
