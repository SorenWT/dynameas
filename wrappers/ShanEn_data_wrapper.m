function [ShanEnOut] = ShanEn_EEG_handle(EEG)

ShanEnOut = zeros(1,EEG.nbchan);

disp(' ')
disp('Computing Shannon Entropy...')

for c = 1:EEG.nbchan
    fprintf([num2str(c) ' ']);
    ShanEnOut(c) = wentropy(EEG.data(c,:),'shannon')/length(EEG.data);
end