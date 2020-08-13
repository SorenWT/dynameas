function [IAFOut] = IAF_EEG_handle_IRASA(spec)

IAFOut = zeros(1,size(spec.osci,2));

disp(' ')
ft_progress('init','text','Computing individual alpha frequency...')

for c = 1:size(spec.osci,2)
        ft_progress(c/dat.nbchan,'Processing channel %d out of %d',c,dat.nbchan);
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
ft_progress('close')