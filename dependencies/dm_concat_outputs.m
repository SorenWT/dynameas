function [newoutputs] = ft_concat_outputs(varargin)

newoutputs = varargin{1};
newoutputs.cfg = {newoutputs.cfg};

for c = 2:length(varargin)
   newoutputs.data = cat(3,newoutputs.data,varargin{c}.data);
   newoutputs.meas = cat(2,newoutputs.meas,varargin{c}.meas);
   newoutputs.cfg{c} = varargin{c}.cfg;
end