function [MFOut] = MF_IRASA_wrapper(EEG,frange,oscifrac)
% Computes the median frequency of the data


if nargin < 2
    frange = [0.5 50];
end


MFOut = zeros(1,EEG.nbchan);

disp(' ')
ft_progress('init','text','Computing median frequency...')

for c = 1:size(EEG.(oscifrac),2)
        ft_progress(c/dat.nbchan,'Processing channel %d out of %d',c,dat.nbchan);
    MFOut(c) = medfreq(EEG.(oscifrac),EEG.freq,frange);
end
ft_progress('close')
