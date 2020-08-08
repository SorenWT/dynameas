%% Separating Fractal and Oscillatory Components in the Power Spectrum of Neurophysiological Signals
% This m-file shows the examples of speparating fractal and oscillatory components from mixed time series. 
%
% Version
%   0.01
%
% Refrerence
%   -Wen H. and Liu Z. (2015) Separating Fractal and Oscillatory Components in 
%    the Power Spectrum of Neurophysiological Signals
%  

%% History
% 0.01 - HGWEN - 12/20/2013 - original file


%% Example of simulation 
% set parameters
srate = 1000;
wlen = 8500;
tlen = 32768;
beta = 1.5;
cons = 20;
Ndata = 2^floor(log2(wlen*0.9));
thetarange = [0 2*pi];
frange = [1 500];

% generate random phase oscilations with specific frequencies or random frequencies
t = 1 : tlen;
f = [10; 20; 40; 60; 80]; % set oscillatory frequency
% f = floor(rand(1,50)*150); % 50 random frequencies
n = length(f); % number of oscilations
theta = rand(1,n)*2*pi;
theta = theta(:);
x = sum(20*cos(2*pi*f*t/srate + repmat(theta,1,tlen)),1);
x = x(:);

% generate fractal time series
sig = amri_sig_genfrac(tlen,Ndata,'beta',beta,'cons',cons,'theta',thetarange);

% add oscillatory components to fractal component
start = ceil(rand(1,1)*(tlen-wlen));
y = sig + x;

% select a time window
Da = y(start:start+wlen-1);

% Separate fractal and oscillatory components
spec = amri_sig_fractal(Da,srate,'frange',frange);

% Display the spectra in log-log scale
figure;
loglog(spec.freq,spec.mixd,'b'); hold on;
loglog(spec.freq,spec.frac,'r'); hold off;
legend('mixd','frac');

%% Example of ECoG data
% set parameter
srate = 1000; % sampling frequency
movingwin = [3 1]; % [window size, sliding step]
frange = [1 100];
win = movingwin(1)*srate;
step = movingwin(2)*srate;

% load ECoG data from one sensor recorded in the left occipital of one
% macaque in eyes-closed awake state, totally 5 mins
load('ECoG_data.mat');

% separate fractal and oscillatory components using sliding window
nwin = floor((length(data) - win)/step);
sig = zeros(win,nwin);
for i = 1 : nwin
    sig(:,i) = data(ceil((i-1)*step)+1 : ceil((i-1)*step)+win);
end
tic
Frac = amri_sig_fractal(sig,srate,'detrend',1,'frange',frange);
Frac.time = (0:step/srate:step*(nwin-1)/srate)';
toc

% fitting power-law function to the fractal power spectra
Frange = [15, 100]; % define frequency range for power-law fitting
Frac = amri_sig_plawfit(Frac,Frange);

% show averaged fractal and oscillatory power spectrum
figure;
subplot(2,1,1);
loglog(Frac.freq,mean(Frac.mixd,2),'b'); hold on
loglog(Frac.freq,mean(Frac.frac,2),'r');
subplot(2,1,2);
plot(Frac.freq, mean(Frac.osci,2));

% show dynamic fractal power spectrum
% figure;
% for i = 1 : length(Frac.time)
%     loglog(Frac.freq,Frac.frac(:,i),'k'); hold on
%     loglog(Frac.Freq, Frac.Plaw(:,i),'r');
%     pause(1);
%     clf;
% end


