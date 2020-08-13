function [SDout] = SD_EEG_handle(EEG)

SDout = [];

ft_progress('init','text','Computing SD...')
for c = 1:EEG.nbchan
        ft_progress(c/dat.nbchan,'Processing channel %d out of %d',c,dat.nbchan);
    SDout(c) = std(EEG.data(c,:));
end
ft_progress('close')