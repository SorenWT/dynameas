function [OsciFrac] = OsciFrac_EEG_wrapper(spec,frange)

OsciFrac = zeros(1,size(spec.osci,2));

if nargin < 2
    frange = [0.5 50];
end

disp(' ')
disp('Computing oscillatory/fractal ratio...')

freqindex = intersect(find(spec.freq > frange(1)),find(spec.freq < frange(2)));
freqs = spec.freq(freqindex);
    
for c = 1:size(spec.osci,2)
    fprintf([num2str(c) ' ']);
    %osci = sgolayfilt(spec.osci(:,c),5,1501);
    osci = spec.osci(freqindex,c);
    %frac = sgolayfilt(spec.frac(:,c),5,1501);
    frac = spec.frac(freqindex,c);
    
    OsciFrac(c) = trapz(freqs,osci)/trapz(freqs,frac);
end