function [asymout] = Frontalasym_data_wrapper(dat,felecs,alphafun)
% This function computes frontal alpha asymmetry using a user-specified
% method
%
% Required inputs: 
%   dat: the raw data, in EEGLAB format
%
% Recommended inputs: 
%   felecs: the electrodes used to calculate asymmetry, given as a cell
%      array of channel labels on the left and right hemispheres (default = {'F4','F3'})
%
% Optional inputs: 
%   alphafun: a function handle to the function and options you want to use
%      to calculate alpha power. (default =
%      @(EEG)Bandpower_data_wrapper(EEG,[8 13]))


if ~exist('felecs','var')
   norm_bandpass = {'F4','F3'}; 
end

if ~exist('alphafun','var')
   alphafun = @(EEG)Bandpower_data_wrapper(EEG,[8 13])
end

alphadat = alphafun(dat);

channames = {dat.chanlocs(:).labels};

for i = 1:2
    if ischar(felecs{i})
        felecs{i} = {felecs{i}};
    end       
    if iscell(felecs{i})
        for ii = 1:length(felecs{i})
            felecs_indx{i}(ii) = find(strcmpi(channames,felecs{i}{ii}));
        end
    elseif isnumeric(felecs{i})
        felecs_indx{i} = felecs{i};
    end
end

asymout = log(alphadat(felecs_indx{1}))-log(alphadat(felecs_indx{2}));

asymout = repmat(asymout,1,dat.nbchan);

