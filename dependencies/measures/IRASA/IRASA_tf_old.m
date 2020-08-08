function [ersp,itc,tpoints,freqs,specs] = IRASA_tf_old(oscifrac,data,srate,frame,varargin)
%the first dimension of data should be times, and the second dimension
%should be trials

% if ~CheckInput(varargin,'TPoints')
%     tpoints = 300;
% end

winsize = 2048*srate/1000 +50;

overlap = 0.98;

[specs,tpoints] = SlidingWindow(@(window)amri_sig_fractal_tf(window,srate),data,winsize,overlap,'mean','off');

outputfreqs = specs(1,1).freq;

% if ~CheckInput(varargin,'frange')
% freqs = [1:50];
% end

for c = 1:length(freqs)
   freqindex(c) = find(abs(outputfreqs-freqs(c)) == min(abs(outputfreqs-freqs(c)))); 
end

%tpoints = linspace(frame(1),frame(2),length(tpoints));
tpoints = frame(1) + tpoints.*1000/srate;

%specs([1:3 (end-2):end]) = [];
%tpoints([1:3 (end-2):end]) = [];

tfdata = [];
tmp = [];
for c = 1:length(specs)
    tfdata = cat(3,tfdata,specs(c).(oscifrac)(freqindex,:));
end

tfdata = permute(tfdata,[1 3 2]);

baseline = find(tpoints < 0);

ersp = mean(bsxfun(@minus,10*log10(abs(tfdata).^2),mean(10*log10(abs(tfdata(:,baseline,:)).^2),2)),3); 
itc = abs(mean((tfdata./abs(tfdata)),3));



















% freqs = intersect(find(specs(1).freq > frange(1)),find(specs(1).freq < frange(2)));
% for c = 1:length(specs)
%     specdata(:,c) = specs(c).(oscifrac)(freqs);
% end
% figure
% mesh(1:length(specs),specs(1).freq(freqs),specdata)
% xlabel('Time')
% ylabel('Frequency')

    
