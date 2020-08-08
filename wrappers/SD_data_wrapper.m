function [SDout] = SD_EEG_handle(EEG)

SDout = [];

disp('Computing SD...')
for c = 1:EEG.nbchan
    fprintf([num2str(c) ' '])
    SDout(c) = std(EEG.data(c,:));
end