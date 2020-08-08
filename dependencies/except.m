function [indxout] = except(indxin,exception)

indxout = indxin(find(~ismember(indxin,exception)));