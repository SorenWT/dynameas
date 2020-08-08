function [IrasaOut] = PLE_IRASA_EEG_wrapper(spec,frange)

disp(' ')
disp('Computing PLE with IRASA...')

if nargin < 2
   frange = [0.5 50]; 
end

tmp = amri_sig_plawfit(spec,frange);
IrasaOut = tmp.Beta;
IrasaOut = IrasaOut';
