function [output] = cfgparse(cfg,field,check)

if nargin > 2
    if isfield(cfg,field) && strcmpi(cfg.(field),check)
        output = 1;
    else
        output = 0;
    end
else
    if isfield(cfg,field)
        output = 1;
    else
        output = 0;
    end
end