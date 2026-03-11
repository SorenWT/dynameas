function [intout] = aper_intercept_fooof_wrapper(dat,varargin)
% Calculates the power-law exponent
% 
% Required inputs:
%   dat: the outputs from a fooof transform
% 
% Optional inputs: optional inputs for this function are in the form of
%   name-value pairs


if nargin < 2
    argsin = {[]};
else
    argsin = varargin;
end

intout = zeros(1,length(dat));

disp(' ')
ft_progress('init','text','Fetching aperiodic intercept...')

for i = 1:length(dat)
   intout(i) = dat(i).aperiodic_params(1); 
end

ft_progress('close')