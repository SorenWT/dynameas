function [MultFracOut] = MultFrac_EEG_handle(EEG)

MultFracOut = zeros(1,EEG.nbchan);

disp(' ')
disp('Computing the spread of Holder exponents...')

for c = 1:EEG.nbchan
    fprintf([num2str(c) ' ']);
    [~,tmp] = dwtleader(EEG.data(c,:));
    MultFracOut(c) = max(tmp)-min(tmp);
end