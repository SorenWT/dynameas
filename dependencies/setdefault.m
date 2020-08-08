function [cfg] = setdefault(cfg,field,value)

if isstruct(cfg)
    if ~isfield(cfg,field)
        cfg.(field) = value;
    end
elseif iscell(cfg)
     if ~CheckInput(cfg,field)
        cfg = [cfg {field value}];
     end 
end