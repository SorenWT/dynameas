function [out] = vert(in)
%makes sure a vector is a column vector

if (isnumeric(in) | iscell(in)) & any(size(in)==1) & (length(size(in))<3)
    if size(in,2) > 1
        out = in';
    else
        out = in;
    end
else
    out = in;
end