function [varout] = parload(filepath,variablename)

if nargin > 1
    tmp = load(filepath);
    varout = tmp.(variablename);
else
    varout = load(filepath);
end