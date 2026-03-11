function fooofdata = transform_fooof(cfg,data)

cfg = setdefault(cfg,'load','yes');
cfg = setdefault(cfg,'save','no');
cfg = setdefault(cfg,'frange',[1 120]);
cfg = setdefault(cfg,'fooofsettings',struct());

disp(' ')
disp('Fitting oscillations and 1/f slope...')
fname = data.filename;


if exist(fullfile(cfg.fldr,[fname '_fooof.mat']),'file') && strcmpi(cfg.load,'yes')
    fooofdata = parload(fullfile(cfg.fldr,[fname '_fooof.mat']),'fooofdata');
else
    %ftdata = eeglab2fieldtrip(data,'preprocessing','none');
    if isfield(data,'nbchan')
        nfft = 2^nextpow2((3/cfg.frange(1))*data.srate);
        for i = 1:data.nbchan
            [pxx(i,:),f] = pwelch(data.data(i,:),[],[],nfft,data.srate);
        end
        
        tic
        for i = 1:data.nbchan
            %try
            %    fooofdata(i) = fooof(f',pxx(i,:),cfg.frange,cfg.fooofsettings);
            %    lnoise{i} = 'yes';
            %catch
                % do all excluding line noise region for now
                nolnoise = except(1:length(f),intersect(find(f>58),find(f<62))); 
                fooofdata(i) = fooof(f(nolnoise)',pxx(i,nolnoise),cfg.frange,cfg.fooofsettings);
                lnoise{i} = 'no';
            %end
        end
        
        for i = 1:data.nbchan
           fooofdata(i).lnoise = lnoise{i};
        end
        toc
        
        for i = 1:length(fooofdata)
            fooofdata(i).filename = data.filename;
            fooofdata(i).f = f;
            fooofdata(i).pxx = pxx(i,:);
            if length(fooofdata(i).aperiodic_params==2)
                fooofdata(i).aper_pxx = 10.^(fooofdata(i).aperiodic_params(1)-fooofdata(i).aperiodic_params(2)*log10(f));
            else % figure this out for the 3-parameter fit
                
            end
            fooofdata(i).resid_pxx = fooofdata(i).pxx-fooofdata(i).aper_pxx;
        end
    else % version for if I do a PSD transform
        
    end
end

if cfgcheck(cfg,'save','yes')
    save(fullfile(cfg.fldr,[fooofdata(1).filename '_fooof.mat']),'fooofdata','-v7.3')
end

end