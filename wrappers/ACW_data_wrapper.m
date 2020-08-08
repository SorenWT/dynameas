function [ACWOut] = ACW_data_wrapper(dat,winsize)
% Calculates the autocorrelation window according to the methods of Honey
% et al. (2012)
% Required inputs: 
%   dat is the raw data, in EEGLAB format
%
% Optional inputs: 
%   winsize is the size of the window in which to compute the ACW (default
%       = 20)
%
% Honey, C. J., Thesen, T., Donner, T. H., Silbert, L. J., Carlson, C. E., 
% Devinsky, O., Doyle, W. K., Rubin, N., Heeger, D. J., & Hasson, U. 
% (2012). Slow Cortical Dynamics and the Accumulation of Information over 
% Long Timescales. Neuron, 76(2), 423?434.


if ~exist('winsize','var')
    winsize = 20;
end

ACWOut = zeros(1,dat.nbchan);

disp(' ')
disp('Computing autocorrelation window...')

for c = 1:dat.nbchan
    fprintf([num2str(c) ' ']);
    ACWOut(c) = ACW_estimation(dat.data(c,:),dat.srate,winsize);
end