function [SpecEnOut] = SpecEn_data_wrapper(dat,frange)
% Computes the spectral entropy of the data
%
% Required inputs: 
%   dat: the raw data in EEGLAB format
% 
% Recommended inputs: 
%   frange: the frequency limits for the spectral entropy - this should be
%       equal to your bandpass at most (default = [1 50])

if ~exist('frange','var')
   frange = [1 50]; 
end

SpecEnOut = zeros(1,dat.nbchan);

disp(' ')
disp('Computing spectral entropy...')

for c = 1:dat.nbchan
    fprintf([num2str(c) ' ']);
    SpecEnOut(c) = pentropy(dat.data(c,:),dat.srate,'Instantaneous',false,'FrequencyLimits',frange);
end