function [MultFracOut] = MultFrac_data_wrapper(EEG)

MultFracOut = zeros(1,EEG.nbchan);

disp(' ')
ft_progress('init','text','Computing the spread of Holder exponents...')

for c = 1:EEG.nbchan
            ft_progress(c/dat.nbchan,'Processing channel %d out of %d',c,dat.nbchan);
    [~,tmp] = dwtleader(EEG.data(c,:));
    MultFracOut(c) = max(tmp)-min(tmp);
end
ft_progress('close')