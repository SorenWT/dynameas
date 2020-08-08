function [specout,specs] = IRASA_window(datain,srate,varargin)
%datain should be in the form channels x time points


if ~CheckInput(varargin,'winsize')
    winsize = round(10*srate);
else
    winsize = round(EasyParse(varargin,'winsize')*srate);
end

if ~CheckInput(varargin,'overlap')
    overlap = 0;
else
    overlap = EasyParse(varargin,'overlap');
end

if ~CheckInput(varargin,'hset')
    hset = [1.1:0.05:1.95 2.05:0.05:2.9];
else
    hset = EasyParse(varargin,'hset');
end

for c = 1:size(datain,1)
   specs(c) = amri_sig_fractal(WindowToMatrix(datain(c,:),winsize,overlap),srate,'hset',hset);
end

% for c = 1:length(specs)
%    specs(c).mixd = nan_mean(specs(c).mixd,2);
%    specs(c).frac = nan_mean(specs(c).frac,2);
%    specs(c).osci = nan_mean(specs(c).osci,2);
% end

specout = specs(1);
specout.mixd = [];
specout.frac = [];
specout.osci = [];
for c = 1:length(specs)
   specout.mixd = cat(2,specout.mixd,nanmean(specs(c).mixd,2));
   specout.frac = cat(2,specout.frac,nanmean(specs(c).frac,2));
   specout.osci = cat(2,specout.osci,nanmean(specs(c).osci,2));
end
end

