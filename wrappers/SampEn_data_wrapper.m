function [SampEnOut] = SampEn_EEG_handle(EEG)

sampEnOut = zeros(1,EEG.nbchan);

disp(' ')
disp('Computing Sample Entropy...')

for c = 1:EEG.nbchan
    fprintf([num2str(c) ' ']);
    tmp = sampen(EEG.data(c,:),2,0.2,1,0);
    sampEnOut(c) = tmp(2);
end