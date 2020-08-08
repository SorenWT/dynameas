function [indx] = Subset_index(cell_all,cell_sub)

if ~iscell(cell_all) && ~iscell(cell_sub)
    cell_all = cellstr(num2str(vert(cell_all)));
    cell_sub = cellstr(num2str(vert(cell_sub)));
end

indx = zeros(1,length(cell_sub));

for c = 1:length(cell_sub)
   indx(c) = find(strcmpi(cell_sub{c},cell_all)); 
end