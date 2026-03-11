function [bp] = Bandpower_data_wrapper(dat,frange,norm_bandpass)
% Computes the power in a given frequency band
%
% Required inputs: 
%   dat: the original data, in datLAB format
%   frange: the frequency range in which you want to compute power, as a
%       tuple
% 
% Recommended inputs:
%   norm_bandpass: determines whether to normalize the power by the total
%       power within some bandwidth. Input a tuple to specify this range,
%       input 'no' to use absolute power (default = 'no')


if ~exist('norm_bandpass','var')
   norm_bandpass = 'no'; 
end

disp(' ')

ft_progress('init','text','Computing power...')
for c = 1:dat.nbchan
    ft_progress(c/dat.nbchan,'Processing channel %d out of %d',c,dat.nbchan);

   [pxx,f] = pwelch(dat.data(c,:),[],[],2^nextpow2((3/2)*dat.srate),dat.srate);
   findx = intersect(find(f > frange(1)),find(f < frange(2)));
   bp(c) = (norm(pxx(findx))^2)./numel(findx);
    %bp(c) = bandpower(dat.data(c,:),dat.srate,frange); 
   if ~strcmpi(norm_bandpass,'no')
       allfindx = intersect(find(f > norm_bandpass(1)),find(f < norm_bandpass(2)));
       bp(c) = bp(c)/((norm(pxx(allfindx))^2)./numel(allfindx));
        %bp(c) = bp(c)/bandpower(dat.data(c,:),dat.srate,norm_bandpass); 
   end
end
ft_progress('close');