function [SampEnOut] = SampEn_EEG_handle(EEG)

sampEnOut = zeros(1,EEG.nbchan);

disp(' ')
ft_progress('init','text','Computing Sample Entropy...')

for c = 1:EEG.nbchan
    ft_progress(c/dat.nbchan,'Processing channel %d out of %d',c,dat.nbchan);
    tmp = sampen(EEG.data(c,:),2,0.2,1,0);
    sampEnOut(c) = tmp(2);
end
ft_progress('close')