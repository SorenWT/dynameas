function [torf] = CheckInput(cellin,name)

torf = 0;

for c = 1:length(cellin)
    if ischar(cellin{c}) && strcmpi(cellin{c},name)
        torf = 1;
        break;
    end
end