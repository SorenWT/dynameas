function [ShanEnOut] = ShanEn_data_wrapper(EEG)

ShanEnOut = zeros(1,EEG.nbchan);

disp(' ')
ft_progress('init','text','Computing Shannon Entropy...')

for c = 1:EEG.nbchan
        ft_progress(c/dat.nbchan,'Processing channel %d out of %d',c,dat.nbchan);
    ShanEnOut(c) = wentropy(EEG.data(c,:),'shannon')/length(EEG.data);
end
ft_progress('close')