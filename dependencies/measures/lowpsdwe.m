function [ple,pdata,freq] = lowpsdwe(time_series,Fs,low_range,high_range,varargin)
% JF_power_law is a modification of Jianfeng Zhang's JF_power_law_nfft1024
% script for calculating the power-law exponent based on Welch's power
% spectral density estimate. The modifications include automatically
% setting the Welch window based on the lowest frequency in the data, and
% adding the option to plot the power spectrum and fit
% 
% Input arguments: 
%    time_series: a vector containing the time series of interest
%    fs: the sampling frequency of the data
%    low_range: the lowest frequency to be included in the fit
%    high_range: the highest frequency to be included in the fit
% 
% Optional key-value pairs: 
%    'Plot','on': plots the power spectrum and the power-law fit
%
% Outputs: 
%    ple: the power-law exponent
%    pdata: the power spectrum
%    freq: the frequencies from the power spectrum estimation

% Parabolic windowing
pwindow = 1-((2.*(1:length(time_series))./(length(time_series)+1)-1).^2);
time_series = pwindow.*time_series;

%endmatching
time_series = time_series-linspace(time_series(1),time_series(2),length(time_series));

%demean
time_series = time_series - mean(time_series);

%from here on same as JF_power_law
nfft = 2^nextpow2((3/low_range)*Fs);

if nfft > length(time_series) || nfft/(2*Fs) < 1
   nfft = []; % if nfft is too large, just use the default Welch window
end

if CheckInput(varargin,'pmethod') && EasyParse(varargin,'pmethod','periodogram')
    [pdata,freq] = periodogram(time_series,[],nfft,Fs); %want 3 cycles of lowest frequency in window
else
    [pdata,freq] = pwelch(time_series,[],[],nfft,Fs); %want 3 cycles of lowest frequency in window
end%     power_spec = psd(HS,time_series,'NFFT',nfft,'Fs',Fs);
%pdata = pxx;
%freq = f;

slope_index = find(freq > low_range & freq < high_range);
%freq = freq(slope_index)';
if CheckInput(varargin,'interp') && EasyParse(varargin,'interp','off')
    linfreq = log10(freq(slope_index));
    fitdata = log10(pdata(slope_index));
else
    fitfreq = log10(freq(slope_index));
    fitdata = log10(pdata(slope_index));
    linfreq = linspace(min(fitfreq),max(fitfreq),length(fitfreq));
    fitdata = interp1(fitfreq,fitdata,linfreq);
end

p = polyfit(vert(linfreq),vert(fitdata),1);
ple = -p(1);


if CheckInput(varargin,'Plot') && EasyParse(varargin,'Plot','on')
    y = p(2) + p(1)*log10(freq(slope_index));
    loglog(freq(slope_index),pdata(slope_index));
    hold on;
    loglog(freq(slope_index),10.^y,'r--');
    xlabel('Log Frequency')
    ylabel('Log Power')
    title(['Estimated PLE is ' num2str(ple)])
end


%[b,bint,r,rint,stats] = regress(log(power_data(slope_index)),[ones(length(power_freq(slope_index)),1) log(power_freq')*p(1)+ p(2)]);
%r_resi = stats(4);
