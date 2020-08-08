function [windowmat,tindx] = WindowToMatrix(vectin,winsize,overlap)

windowmat = [];
indx = 1;

if size(vectin,1) == 1
    vectin = vectin';
end

tindx = [];
tindx = [tindx indx];

while (indx+winsize) < length(vectin)
    windowmat = cat(2,windowmat,vectin(indx:(indx+winsize-1)));
    indx = indx+winsize*(1-overlap);
    tindx = [tindx indx];
end