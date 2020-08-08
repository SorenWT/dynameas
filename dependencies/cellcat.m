function [cellout] = cellcat(stringin,cellin,joiner,right)

if nargin < 4
    right = 0;
end

if right == 0
    for c = 1:length(cellin)
        if iscell(cellin{c})
            cellout{c} = cellcat(stringin,cellin{c},joiner,0);
        else
            cellout{c} = [stringin joiner cellin{c}];
        end
    end
elseif right == 1
    for c = 1:length(cellin)
        if iscell(cellin{c})
            cellout{c} = cellcat(stringin,cellin{c},joiner,1);
        else
            cellout{c} = [cellin{c} joiner stringin];
        end
    end
end
end