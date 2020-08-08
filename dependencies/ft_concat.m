function [data] = ft_concat(data)

tmp = cat(2,data.trial{:});
data.trial = [];
data.trial{1} = tmp;
data.time = [];
data.time{1} = linspace(0,size(data.trial{1},2)/data.fsample,size(data.trial{1},2));