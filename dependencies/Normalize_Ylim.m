function [newlim] = Normalize_Ylim(fighandle,equal)

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
        Ylims(c,:) = get(ax(c),'YLim');
    catch
    end
end

newlim = [min(Ylims(:,1)) max(Ylims(:,2))];

if equal
    maxabs = max(abs(newlim));
    newlim = [-maxabs maxabs];
end
for c = 1:length(ax)
    try
        set(ax,'YLim',newlim)
    catch
    end
end