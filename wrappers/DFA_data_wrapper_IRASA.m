function [dfaout] = DFA_EEG_wrapper_IRASA(spec,oscifrac,frange)

dfaout = zeros(1,length(spec));

if nargin < 2
    oscifrac = 'osci';
end

if nargin < 3
    frange = [8 13];
end

disp(' ')
disp('Computing DFA exponent...')

for c = 1:length(dfaout)
    fprintf([num2str(c) ' '])
    findex = intersect(find(spec(c).freq > frange(1)),find(spec(c).freq < frange(2)));
    ts = sqrt(trapz(spec(c).freq(findex),spec(c).(oscifrac)(findex,:),1));
    [~,~,tmp] = FMF(ts,nextpow2(spec(1).srate),nextpow2(0.5*length(ts))/2,50,2,1,'Plot','off');
    dfaout(c) = tmp.MF;
end
