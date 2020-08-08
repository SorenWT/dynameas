function [miout] = MI_data_wrapper(data,lf,hf,varargin)

siz = size(lf);
if siz(1) == 1 || siz(2) == 1
    lf = horz(lf);
end

siz = size(hf);
if siz(1) == 1 || siz(2) == 1
    hf = horz(hf);
end

disp(' ')
disp('Computing modulation index...')

ftdat = eeglab2fieldtrip(data,'preprocessing','none');

%for i = 1:size(lf,1)
cfg = []; cfg.bpfilter = 'yes'; cfg.bpfreq = lf; cfg.hilbert = 'angle'; %cfg.bpfilttype = 'fir';
lfdat = ft_preprocessing(cfg,ftdat);
%lfdat = ft_cont2epoch(lfdat,4);

%end

%for i = 1:size(lf,1)
cfg = []; cfg.bpfilter = 'yes'; cfg.bpfreq = hf; cfg.hilbert = 'abs';
hfdat = ft_preprocessing(cfg,ftdat);
%hfdat = ft_cont2epoch(hfdat,4);
%end

if ~CheckInput(varargin,'nbins')
    nbins = 18;
else
    nbins = EasyParse(varargin,'nbins');
end

for c = 1:data.nbchan
    fprintf([num2str(c) ' ']);
    tmp = get_mi(lfdat.trial{1}(c,:)',hfdat.trial{1}(c,:)',nbins);%,200,3);
    miout(c) = tmp.MI;
    %miout(c) = tmp.MIp;
end