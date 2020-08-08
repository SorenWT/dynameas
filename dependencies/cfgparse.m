function [output] = cfgparse(cfg,field)

if isfield(cfg,field)
    output = cfg.(field);
else
    output = NaN;
end