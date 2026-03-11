function [IrasaOut] = aper_intercept_irasa_wrapper(spec,frange)

disp(' ')
disp('Computing aperiodic intercept with IRASA...')

if nargin < 2
   frange = [0.5 50]; 
end

tmp = amri_sig_plawfit(spec,frange);
IrasaOut = tmp.Cons;
IrasaOut = IrasaOut';
