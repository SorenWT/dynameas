function [IAFOut] = IAF_EEG_handle(EEG)

IAFOut = zeros(1,EEG.nbchan);

disp(' ')
disp('Computing individual alpha frequency...')
%EEG = pop_eegfiltnew(EEG, 7, 13, 826, 0, [], 0);
[spectra,freqs] = spectopo(EEG.data,0,EEG.srate,'plot','off');
for c = 1:EEG.nbchan
    fprintf([num2str(c) ' ']);
    peaks = findpeaks(spectra(c,:));
    if ~isempty(peaks)
        peaks = peaks.loc;
        peaks = peaks(intersect(find(peaks >= 7),find(peaks <= 15)));
        IAFOut(c) = freqs(find(spectra(c,:) == max(spectra(c,peaks))));
    else
        IAFOut(c) = NaN;
    end
end