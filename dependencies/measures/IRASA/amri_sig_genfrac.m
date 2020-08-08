%% Simulate fractal time series generating function
%   amri_sig_genfrac()
%
% Usage
%   frac = amri_sig_genfrac(N,Nfft,...)
% 
% Inputs
%   N    - The length of output fractal time series. It would be slow if N is very large.
%   Nfft - The number of point for fft analysis
%
% Outputs
%   frac - Generated fractal time series, a column verctor.
%
% Keywords
%   theta - [thetamin thetamax], the range of phase uniform distribution
%   beta  - power-law exponent
%   cons  - intercept of 1/f^beta line
% 
% See also
%   amri_sig_fractal
%
% Version
%   0.01
% 
% Reference
%   -Wen H. and Liu Z. Separating Fractal and Oscillatory Components in 
%    the Power Spectrum of Neurophysiological Signals

%% history
% 0.01 - HGWEN - 12/20/2013 - original file


%%
function frac = amri_sig_genfrac(N,Nfft,varargin)
if nargin<2
    eval('help amri_sig_genfrac');
    return
end

%% Defaults
thetarange = [0 2*pi];
beta = 1;
cons = 1;

%% Keywords
for i = 1:2:size(varargin,2) 
    Keyword = varargin{i};
    Value   = varargin{i+1};
    if strcmpi(Keyword,'theta')
        thetarange = Value;
    elseif strcmpi(Keyword,'beta')
        beta = Value;
    elseif strcmpi(Keyword,'cons')
        cons = Value;
    else
        warning(['amri_sig_genfrac(): unknown keyword ' Keyword]);
    end
end

%% generate fractal time series
t = 1 : N;
k = (1 : Nfft/2)';
kt = k*t;
theta = unifrnd(thetarange(1),thetarange(2),Nfft/2,1);
Ck = repmat(sqrt(cons*k.^(-beta)*Nfft),1,N(1));
frac = sum(Ck.*cos(2*pi*kt/Nfft - repmat(theta,1,N)));
frac = frac(:);
end