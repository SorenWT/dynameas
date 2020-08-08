function [out] = horz(in)
%makes sure a vector is a row vector

if size(in,1) > 1
   out = in';
else
    out = in;
end