function [out] = vert(in)
%makes sure a vector is a column vector

if size(in,2) > 1
   out = in';
else
    out = in;
end