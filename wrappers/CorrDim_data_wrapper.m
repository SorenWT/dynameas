function [CorrDimOut] = CorrDim_EEG_handle(EEG)

CorrDimOut = zeros(1,EEG.nbchan);

disp(' ')
disp('Computing correlation dimension...')

for c = 1:EEG.nbchan
    fprintf([num2str(c) ' ']);
    CorrDimOut(c) = correlationDimension(EEG.data(c,:));
end