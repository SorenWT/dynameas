function manlegend(labels,colors)

g = gca; 
allobj = findall(g.Children,'HandleVisibility','on');

for c = 1:length(allobj)
   set(allobj(c),'HandleVisibility','off') 
end

if ~iscell(colors)
    for c = 1:size(colors,1)
        patch([NaN NaN],[NaN NaN],colors(c,:))
    end
else
    for c = 1:length(colors)
        patch([NaN NaN],[NaN NaN],colors{c})
    end
end

legend(labels)