function c = palecol(c,palefact)

if nargin < 2
   palefact = 0.5; 
end

t = [1 1 1];
d = t - c;
c = c + (d * palefact);