function [fwindows,tpoints,meanwindow] = SlidingWindow(funchandle,data,winsize,overlap,varargin)
%applies a function in a sliding window. If input is a matrix, the sliding
%window will be applied columnwise - the function handle must also operate
%across the first dimension, otherwise results will be incorrect

if size(data,1) == 1
   data = data'; 
end

ii=1;   % Windows counter
while true
    % Begining and ending of the current window
    SWindow=[round(1+(ii-1)*winsize*(1-overlap)), round(winsize+(ii-1)*winsize*(1-overlap))];
    % Check if index exceeds vector dimensions. If so, break!
    if SWindow(2)>=size(data,1)
        break
    end
    % ACF computation into the window (normalized between -1 and 1)
    
    fwindows(ii,:) = funchandle(data(SWindow(1):SWindow(2),:));
    %tpoints(ii) = SWindow(1);
    tpoints(ii) = round(mean(SWindow));
    % Next window
    ii=ii+1;
end

if ~EasyParse(varargin,'mean','off')
    meanwindow = mean(fwindows);
end