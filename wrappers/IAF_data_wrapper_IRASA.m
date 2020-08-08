function [IAFOut] = IAF_EEG_handle_IRASA(spec)

IAFOut = zeros(1,size(spec.osci,2));

disp(' ')
disp('Computing individual alpha frequency...')

for c = 1:size(spec.osci,2)
    fprintf([num2str(c) ' ']);
    osci = spec.osci(:,c);
    %osci = sgolayfilt(spec.osci(:,c),5,1501);
    [~,peaks] = findpeaks(osci);
    if ~isempty(peaks)
        %peaks = peaks.loc;
        peaks = peaks(intersect(find(spec.freq(peaks) >= 7),find(spec.freq(peaks) <= 15)));
        IAFOut(c) = spec.freq(find(osci == max(osci(peaks))));
    else
        IAFOut(c) = NaN;
    end
end