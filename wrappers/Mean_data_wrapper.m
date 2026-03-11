function mdat = Mean_data_wrapper(dat)
% Computes the mean (not useful for EEG, more for source-level data or EMG)

% Required inputs:
%   dat: raw data, in EEGLAB format

for i = 1:dat.nbchan
    mdat(i) = nanmean(dat.data(i,:));
end