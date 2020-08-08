function [dfaout] = DFA_EEG_wrapper(EEG)

dfaout = zeros(1,EEG.nbchan);

disp(' ')
disp('Computing DFA exponent...')

for c = 1:EEG.nbchan
   [~,~,tmp] = FMF(EEG.data(c,:),nextpow2(EEG.srate),nextpow2(0.5*size(EEG.data,2))/2,50,2,1);
   dfaout = tmp.MF;
end