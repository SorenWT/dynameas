function [newlim] = Normalize_Clim(fighandle,equal)

if nargin < 2
    equal = 0;
end

if isa(fighandle,'matlab.ui.Figure')
ax = findall(fighandle,'Type','Axes');
elseif isa(fighandle,'matlab.graphics.axis.Axes')
    ax = fighandle;
end

for c = 1:length(ax)
    try
        Clims(c,:) = get(ax(c),'CLim');
    catch
    end
end

newlim = [min(Clims(:,1)) max(Clims(:,2))];

if equal
    maxabs = max(abs(newlim));
    newlim = [-maxabs maxabs];
end
for c = 1:length(ax)
    try
        set(ax,'CLim',newlim)
    catch
    end
end