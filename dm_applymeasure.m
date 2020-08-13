function [outputs] = dm_applymeasure(cfg)
% ft_applymeasure applies a certain measurement (ex. PLE, DFA, LZC) to a
% data set
%
% Input arguments:
%
% cfg: a configuration structure with the following fields and defaults
%
%      Basic options:
%
%      measure: a cell array containing the function handles to apply to
%         the data (required input). These function handles should take an
%         EEG structure, or at least a structure with fields 'data' (in the
%         format channels x samples) and 'srate' (the sampling rate).
%      format: 'eeglab' or 'fieldtrip' (default = 'eeglab');
%      dir: the base directory (default = get from user interface)
%      files: a string telling what files to select. Usually this
%         includes a wildcard character ('*'). For example, to select all
%         files in the folder ending in 'myfiles.mat', enter '*myfiles.mat'.
%         You can also set up nested directory structures this way. To
%         select all files named 'myfiles.mat' in all directories one level
%         down from the current one, enter */*myfiles.mat. (default =
%         '*.mat' for fieldtrip format, '*.set' for eeglab format)
%      outfile: the file into which to write the output. Can specify 'none'
%         to skip writing to disk (default = 'outputs.mat')
%
%
%      transform: a cell array of configuration structures for transforming
%         the data prior to applying measures. Each config structure should
%         contain the following fields:
%             func: the transform function to apply, as a function handle
%             fldr: the folder in which to load or save the transformed
%               data, if applicable (defaults to the same folder as the
%               original file
%         See each transform function's documentation to see the config
%         options specific to that function
%
%      Because of the comptutationally expensive nature of methods like the
%      IRASA, data transformations apply to the whole set of measures.
%      Currently, you cannot apply a filter or amplitude envelope with one
%      measure and not with others. To facilitate batch processing given
%      this constraint, you can input a cell array of cfg structures
%      instead of a single structure, each with their own set of
%      transformation options. The final outputs will be concatenated
%      together.
%
%      Other options:
%
%      continue: This option loads the output file and restarts the loop
%         from the last subject saved. Cannot be used with a cell array cfg
%         input (default = 'no')
%      ftvar: for fieltrip format - the variable name in the .mat file
%         which corresponds to your data (default = the first structure
%         found in the data file)
%      concatenate: concatenate resting-states organized into trials
%         (default = 'yes')
%      parallel: a structure with the following options:
%         do_parallel: use MATLAB parallel processing (default = 'no');
%         pool: size of parallel pool (default = system default pool size)
%      single: convert fieldtrip data to single precision to save memory
%         (default = 'no')
%      subsrange: calculate on only a given range of files - for use with
%         clusters like ComputeCanada (default = 1:length(files))
%      readonly: option in case the relevant files are in a read-only
%         location - saves intermediate outputs like IRASA spectra or
%         envelope data in the directory of the output file, not the input
%         files (default = 'no')
%
%
% Outputs:
%
% outputs: a structure with the following fields:
%      data: the applied measures
%      dimord: a string telling you the order of the dimensions of 'data'.
%         For example, 'sub_chan_meas' means that the first dimension is
%         subjects, the second dimension is channels, and the third
%         dimension is measures.
%      sub: the file names that the function read
%      chan: the channel labels of the data files
%      meas: the measures you input in cfg.measure
%      elec or grad: the electrode or gradiometer position information -
%         used later in ft_measurestatistics
%      chanlocs: if eeglab data were input, this is the chanlocs structure
%         of the data - to be used later in ft_measurestatplot
%      cfg: the cfg structure used by the function, including the various
%         defaults that were set


if iscell(cfg)
    alloutputs = cell(size(cfg));
    for i = 1:length(cfg)
        alloutputs{i} = dm_applymeasure(cfg{i});
    end
    outputs = dm_concat_outputs(alloutputs{:});
    outputs.cfg = cfg;
    parsave(cfg{1}.outfile,'outputs',outputs);
else
    %% Set up defaults
    cfg = setdefault(cfg,'format','eeglab');
    
    if ~cfgcheck(cfg,'files')
        if cfgcheck(cfg,'format','eeglab')
            cfg.files = '*.set';
        else
            cfg.files = '*.mat';
        end
    end
    
    cfg = setdefault(cfg,'writeouts','yes');
    
    cfg = setdefault(cfg,'outfile',fullfile(pwd,'outputs.mat'));
    if cfg.outfile(1) ~= filesep && cfg.outfile(2) ~= ':'  && strcmpi(cfg.writeouts,'yes') %check for relative path
        cfg.outfile = fullfile(pwd,cfg.outfile);
    end
    
    cfg = setdefault(cfg,'transform',{});
    
    if ~isempty(cfg.transform)
        tmp = tokenize(cfg.outfile,filesep);
        fldr = fullfile(tmp{1:end-1});
        if fldr(1) ~= filesep
            fldr = [filesep fldr];
        end
        for i = 1:length(cfg.transform)
            cfg.transform{i}.fldr = fldr;
        end
    end
    
    cfg = setdefault(cfg,'continue','no');
    
    cfg = setdefault(cfg,'concatenate','yes');
    
    cfg = setdefault(cfg,'single','no');
    
    if any(strcmpi(cfg.format,{'afni','nifti','spm','analyze'}))
        cfg = setdefault(cfg,'writeimage','yes');
    else
        cfg = setdefault(cfg,'writeimage','no');
    end
    
    if ~cfgcheck(cfg,'parallel')
        cfg.parallel.do_parallel = 'no';
    end
    
    cfg.parallel = setdefault(cfg.parallel,'pool','default');
    
    %% Start up fieldtrip, eeglab if necessary
    if cfgcheck(cfg,'format','eeglab')
        %        eeglab rebuild
    end
    
    ft_defaults
    
    ftdir = which('ft_defaults');
    ftdir = char(extractBefore(ftdir,'ft_defaults.m'));
    addpath(fullfile(ftdir,'external','eeglab'))
    
    %% Find the base directory
    dirname = cfgparse(cfg,'dir');
    if isnan(dirname)
        dirname = uigetdir;
    end
    
    files = dir(fullfile(dirname,cfg.files));
    
    %% Set up the outputs structure
    if cfgcheck(cfg,'continue','yes')
        load(cfg.outfile)
    else
        outputs = struct;
        outputs.dimord = 'sub_chan_meas';
        
        outputs.meas = cfg.measure;
        outputs.startsub = 1;
        
    end
    
    cfg = setdefault(cfg,'subsrange',outputs.startsub:length(files));
    
    % Trim the sub range to be no larger than the number of files to avoid error
    cfg.subsrange = cfg.subsrange(find(cfg.subsrange <= length(files)));
    
    
    if cfgcheck(cfg.parallel,'do_parallel','no')
        
        for i = cfg.subsrange
            
            %% Load the data
            outputs.sub{i} = files(i).name;
            outputs.startsub = i+1;
            
            disp(' ')
            disp(['Now processing ' files(i).name])
            
            
            switch cfg.format
                case 'fieldtrip'
                    EEG = struct;
                    allvars = load(fullfile(files(i).folder,files(i).name));
                    if ~cfgcheck(cfg,'ftvar')
                        names = fieldnames(allvars);
                        for c = 1:length(names)
                            if isstruct(allvars.(names{c}))
                                data = allvars.(names{c});
                                clear allvars
                                break
                            end
                        end
                        clear names
                    else
                        data = allvars.(ftvar);
                    end
                    
                    if cfgcheck(cfg,'concatenate','yes')
                        data = ft_concat(data);
                    end
                    
                    if i == outputs.startsub-1
                        if isfield(data,'label')
                            outputs.chan = data.label;
                        end
                        if isfield(data,'elec')
                            outputs.elec = data.elec;
                        elseif isfield(data,'grad')
                            outputs.grad = data.grad;
                        end
                    end
                    
                    EEG = ft2eeglab(data);
                case 'eeglab'
                    EEG = pop_loadset( 'files(i).name', files(i).name, 'filepath', files(i).folder);
                    outputs.chanlocs = EEG.chanlocs;
                    if i == outputs.startsub-1
                        data = eeglab2fieldtrip(EEG,'preprocessing','none');
                        
                        if cfgcheck(cfg,'concatenate','yes')
                            data = ft_concat(data);
                            EEG = ft2eeglab(data);
                            EEG.chanlocs = outputs.chanlocs;
                        end
                        
                        if isfield(data,'label')
                            outputs.chan = data.label;
                        end
                        if isfield(data,'elec')
                            outputs.elec = data.elec;
                        elseif isfield(data,'grad')
                            outputs.grad = data.grad;
                        end
                    end
                case 'afni'
                    opt = struct; opt.format = 'vector';
                    [~,brik,brikinfo] = BrikLoad(fullfile(files(i).folder,files(i).name));
                    EEG = eeg_emptyset();
                    EEG.data = brik;
                    brik = [];
                    EEG.srate = 1000/brikinfo.TAXIS_FLOATS(2);
                    EEG.etc.dims = brikinfo.DATASET_DIMENSIONS(1:3);
                    EEG = eeg_checkset(EEG);
                    EEG = ft_struct2single(EEG);
                    
                    if i == outputs.startsub-1
                        outputs.ext{i} = char(extractAfter(files(i).name,'+'));
                        outputs.chan = cellcat('vox',cellstr(num2str([1:EEG.nbchan]')),'',0);
                        outputs.hdr = brikinfo;
                    end
                case 'nifti'
                    info = niftiinfo(fullfile(files(i).folder,files(i).name));
                    
                    data = niftiread(fullfile(files(i).folder,files(i).name));
                    EEG = eeg_emptyset();
                    EEG.data = reshape(data,[],info.ImageSize(4));
                    data = [];
                    if strcmpi(info.TimeUnits,'second')
                        EEG.srate = 1/info.ImageSize(4);
                    elseif strcmpi(info.TimeUnits,'millisecond')
                        EEG.srate = 1000/info.ImageSize(4);
                    else
                        error('Dynameas error: unknown time unit in nifti header')
                    end
                    
                    EEG.etc.dims = info.ImageSize;
                    EEG = eeg_checkset(EEG);
                    EEG = ft_struct2single(EEG);
                    
                    if i == outputs.startsub-1
                        outputs.hdr = info;
                        outputs.chan = cellcat('vox',cellstr(num2str([1:EEG.nbchan]')),'',0);
                    end
                case {'spm','analyze'}
                    mkdir tmp
                    gunzip(fullfile(files(i).folder,files(i).name),fullfile(pwd,tmp))
                    cd tmp
                    tmpfiles = dir('*.img');
                    V = spm_vol(extractfield(tmpfiles,'name'));
                    V = cat(1,V{:});
                    dat = spm_read_vols(V);
                    EEG = eeg_emptyset(EEG);
                    EEG.data = reshape(dat,[],length(V));
                    dat = [];
                    EEG.srate = 1; %fix this later
                    EEG.etc.dims = V(1).dim;
                    EEG = eeg_checkset(EEG);
                    EEG = ft_struct2single(EEG);
                    
                    if i == outputs.startsub-1
                        outputs.hdr = V;
                        outputs.chan = cellcat('vox',cellstr(num2str([1:EEG.nbchan]')),'',0);
                    end
                    cd ..
                    system('rm -r tmp') % fix later - requires linux or mac
                case {'raw','mne'}
                    cfg = []; cfg.dataset = fullfile(files(i).folder,files(i).name);
                    data = ft_preprocessing(cfg,data);
                    
                    if cfgcheck(cfg,'concatenate','yes')
                        data = ft_concat(data);
                    end
                    
                    if i == outputs.startsub-1
                        if isfield(data,'label')
                            outputs.chan = data.label;
                        end
                        if isfield(data,'elec')
                            outputs.elec = data.elec;
                        elseif isfield(data,'grad')
                            outputs.grad = data.grad;
                        end
                    end
                    
                    EEG = ft2eeglab(data);
                otherwise
                    error('Dynameas error: format not recognized. Please check that cfg.format is of a known type')
            end
            
            if i == outputs.startsub-1
                if cfgcheck(cfg,'mask')
                    if contains(cfg.mask,'nii')
                        mask = niftiread(cfg.mask);
                    elseif contains(cfg.mask,'brik')
                        [~,mask,tmp] = BrikLoad(cfg.mask);
                    elseif contains(cfg.mask,'img')
                        tmpv = spm_vol(cfg.mask);
                        mask = spm_read_vols(tmpv);
                    end
                    mask = reshape(logical(mask),[],1);
                elseif cfgcheck(cfg,'channel')
                    mask = cfg.channel;
                end
                outputs.mask = mask;
            end
            EEG = pop_select(EEG,'channel',mask);
            
            
            if cfgcheck(cfg,'single','yes')
                data = ft_struct2single(EEG);
            end
            
            EEG.files(i).name = files(i).name;
            EEG.filepath = files(i).folder;
            
            
            %% Apply the transforms and measures
            for c = 1:length(cfg.transform)
                transfunc = cfg.transform{c}.func;
                EEG = transfunc(cfg.transform{c},EEG);
            end
            
            if strcmpi(cfg.writeouts,'yes')
                outputs.data = zeros(length(cfg.subsrange),length(mask),length(cfg.measure));
                for c = 1:length(cfg.measure)
                    tmpdata = cfg.measure{c}(EEG);
                    outputs.data(i,horz(mask),c) = tmpdata;
                end
            else
                outdata = zeros(1,length(mask),length(cfg.measure));
                outdata(1,horz(mask),c) = tmpdata;
            end
            
            outputs.cfg = cfg;
            
            %% Save after every subject so you can continue later
            try
                if strcmpi(cfg.writeouts,'yes')
                    save(cfg.outfile,'outputs','-v7.3');
                end

                if any(strcmpi({'afni','nifti','spm','analyze'},cfg.format)) && strcmpi(cfg.writeimage,'yes')
                    % write an image containing the outputs
                    switch cfg.format
                        case 'afni'
                            currdir = pwd;
                            outs = tokenize(cfg.outfile,'/');
                            outdir = fullfile(outs{1:end-1});
                            if ~isempty(outdir)
                                cd(outdir)
                            end
                            brikout = reshape(squeeze(outdata(1,:,:)),outputs.hdr.DATASET_DIMENSIONS(1),...
                                outputs.hdr.DATASET_DIMENSIONS(2),outputs.hdr.DATASET_DIMENSIONS(3),length(outputs.meas));
                            brikinfoout = outputs.hdr;
                            brikinfoout.DATASET_DIMENSIONS(4) = length(outputs.meas);
                            tmp = cellfun(@func2str,outputs.meas,'UniformOutput',false);
                            tmp = extractBefore(tmp,'_wrapper'); tmp = erase(tmp,'_');
                            tmp = cellcat('~',tmp,'',1); tmp = cat(2,tmp{:}); tmp(end) = [];
                            brikinfoout.BRICK_LABS = tmp;
                            brikinfoout.BRICK_TYPES = [1]; %store data as shorts
                            brikinfoout.BRICK_STATS = []; %automatically set
                            brikinfoout.BRICK_FLOAT_FACS = [];%automatically set
                            brikinfoout.IDCODE_STRING = '';%automatically set
                            brikinfoout.RootName = ''; %that'll get set by WriteBrik
                            brikinfoout.DATASET_RANK(1) = length(outputs.meas);
                            % Finish fixing this later so the output
                            % makes sense to afni
                            
                            opts.Scale = 1;
                            suf = outs{end}; suf = erase(suf,'.mat');
                            subpref = extractBefore(outputs.sub{i},'+');
                            opts.Prefix = [subpref '_' suf];
                            opts.verbose = 0;
                            WriteBrik(brikout,brikinfoout,opts)
                            cd(currdir)
                        case 'nifti'
                            niftiout = reshape(squeeze(outdata(1,:,:)),outputs.hdr.ImageSize(1),...
                                outputs.hdr.ImageSize(2),outputs.hdr.ImageSize(3),length(outputs.meas));
                            headerout = outputs.hdr;
                            headerout.hdr.ImageSize(4) = length(outputs.meas);
                            headerout.Filename = [char(extractBefore(headerout.Filename,'nii')) '_' outs{end} '.nii'];
                            headerout.Filemoddate = char(datetime);
                            headerout.PixelDimensions(4) = 1;
                            tmp = cellfun(@func2str,outputs.meas,'UniformOutput',false);
                            tmp = extractBefore(tmp,'_wrapper'); tmp = erase(tmp,'_');
                            tmp = cellcat('~',tmp,'',1); tmp = cat(2,tmp{:}); tmp(end) = [];
                            headerout.Description = tmp;
                            niftiwrite(niftiout,headerout.Filename,headerout);
                        case {'spm','analyze'}
                            for q = 1:length(outputs.meas)
                                hdr = outputs.hdr(1);
                                hdr.dt = [16 0];
                                outdata = squeeze(single(outdata(:,:,q)));
                                outdata = reshape(outdata,hdr.dim(1),hdr.dim(2),hdr.dim(3));
                                hdr.fname = [char(extractBefore(hdr.fname,'.img')) '_' outs{end} '_meas' num2str(q) '.img'];
                                spm_write_vol(hdr,outdata);
                            end
                            tar([char(extractBefore(hdr.fname,'_meas')) '.tar.gz'],'*meas*.img')
                            system(['rm *meas*.img']) % again, only works for linux and mac
                    end
                end  
                disp(['Results saved in ' cfg.outfile])
            catch
                warning(['dynameas warning: Saving failed. Try writing to a different directory,' newline ...
                    'or request output argument of dm_applymeasure rather than writing to file'])
            end
        end
        
        disp('dm_applymeasure finished')
        disp('')
    else
        %% Parallel version
        
        if cfgcheck(cfg,'format','fieldtrip')
            allvars = parload(fullfile(files(1).folder,files(1).name));
            if ~cfgcheck(cfg,'ftvar')
                names = fieldnames(allvars);
                for c = 1:length(names)
                    if isstruct(allvars.(names{c}))
                        tmpdata = allvars.(names{c});
                        allvars = [];
                        break
                    end
                end
                names = [];
            else
                tmpdata = allvars.(ftvar);
            end
        else
            EEG = pop_loadset('files(i).name',files(1).name,'filepath',files(1).folder);
            outputs.chanlocs = EEG.chanlocs;
            tmpdata = eeglab2fieldtrip(EEG,'preprocessing','none');
        end
        
        if isfield(tmpdata,'label')
            outputs.chan = tmpdata.label;
        end
        if isfield(tmpdata,'elec')
            outputs.elec = tmpdata.elec;
        elseif isfield(tmpdata,'grad')
            outputs.grad = tmpdata.grad;
        end
        
        currpool = gcp('nocreate');
        if ~cfgcheck(cfg.parallel,'pool','default')
            if ~isempty(currpool) && (currpool.NumWorkers ~= cfg.parallel.pool)
                delete(gcp('nocreate'))
                parpool(cfg.parallel.pool)
            elseif isempty(currpool)
                parpool(cfg.parallel.pool);
            end
        end
        
        clear data
        
        sub = cell(1,length(files));
        outdata = cell(1,length(files));
        
        
        parfor i = 1:length(cfg.subsrange)
            
            %% Load the data
            files(i).name = files(cfg.subsrange(i)).name;
            sub{i} = files(cfg.subsrange(i)).name;
            %outputs.startsub = i+1;
            
            disp(' ')
            disp(['Now processing subject ' num2str(i)])
            
            
            if cfgcheck(cfg,'format','fieldtrip')
                EEG = struct;
                allvars = parload(fullfile(files(cfg.subsrange(i)).folder,files(i).name));
                if ~cfgcheck(cfg,'ftvar')
                    names = fieldnames(allvars);
                    for c = 1:length(names)
                        if isstruct(allvars.(names{c}))
                            data = allvars.(names{c});
                            allvars = [];
                            break
                        end
                    end
                    names = [];
                else
                    data = allvars.(ftvar);
                end
                
                if cfgcheck(cfg,'concatenate','yes')
                    data = ft_concat(data);
                end
                
                EEG = ft2eeglab(data);
            else
                EEG = pop_loadset( 'files(i).name', files(i).name, 'filepath', files(cfg.subsrange(i)).folder);
            end
            
            EEG.files(i).name = files(cfg.subsrange(i)).name;
            EEG.filepath = files(cfg.subsrange(i)).folder;
            
            %% Apply the transforms and measures
            for c = 1:length(cfg.transform)
                transfunc = cfg.transform{c}.func;
                EEG = transfunc(cfg.transform{c},EEG);
            end
            
            
            for c = 1:length(cfg.measure)
                outdata{i}(1,:,c) = cfg.measure{c}(EEG);
            end
            
            
            
        end
        
        outputs.data = cat(1,outdata{:});
        outputs.sub = sub;
        outputs.cfg = cfg;
        
        if ~cfgcheck(cfg,'outfile','none')
            try
                parsave(cfg.outfile,'outputs',outputs,'-v7.3');
                
                if any(strcmpi({'afni','nifti','spm','analyze'},cfg.format)) && strcmpi(cfg.writeimage,'yes')
                    % write an image containing the outputs
                    switch cfg.format
                        case 'afni'
                            currdir = pwd;
                            outs = tokenize(cfg.outfile,'/');
                            outdir = fullfile(outs{1:end-1});
                            if ~isempty(outdir)
                                cd(outdir)
                            end
                            for i = 1:length(outputs.sub)
                                brikout = reshape(squeeze(outputs.data(i,:,:)),outputs.hdr.DATASET_DIMENSIONS(1),...
                                    outputs.hdr.DATASET_DIMENSIONS(2),outputs.hdr.DATASET_DIMENSIONS(3),length(outputs.meas));
                                brikinfoout = outputs.hdr;
                                brikinfoout.DATASET_DIMENSIONS(4) = length(outputs.meas);
                                tmp = cellfun(@func2str,outputs.meas,'UniformOutput',false);
                                tmp = extractBefore(tmp,'_wrapper'); tmp = erase(tmp,'_');
                                tmp = cellcat('~',tmp,'',1); tmp = cat(2,tmp{:}); tmp(end) = [];
                                brikinfoout.BRICK_LABS = tmp;
                                brikinfoout.BRICK_TYPES = [1]; %store data as shorts
                                brikinfoout.BRICK_STATS = []; %automatically set
                                brikinfoout.BRICK_FLOAT_FACS = [];%automatically set
                                brikinfoout.IDCODE_STRING = '';%automatically set
                                brikinfoout.RootName = ''; %that'll get set by WriteBrik
                                brikinfoout.DATASET_RANK(1) = length(outputs.meas);
                                % Finish fixing this later so the output
                                % makes sense to afni
                                
                                opts.Scale = 1;
                                suf = outs{end}; suf = erase(suf,'.mat');
                                subpref = extractBefore(outputs.sub{i},'+');
                                opts.Prefix = [subpref '_' suf];
                                opts.verbose = 0;
                                WriteBrik(brikout,brikinfoout,opts)
                            end
                            cd(currdir)
                        case 'nifti'
                            niftiout = reshape(outputs.data,outputs.hdr.ImageSize(1),...
                                outputs.hdr.ImageSize(2),outputs.hdr.ImageSize(3),length(outputs.meas));
                            
                    end
                end
                
                disp('dm_applymeasure finished')
                disp('')
                disp(['Results saved in ' cfg.outfile])
            catch
                warning(['dynameas warning: Saving failed. Try writing to a different directory,' newline ...
                    'or request output argument of dm_applymeasure rather than writing to file'])
            end
        end
        
    end
end

end

%%
function [EEG] = SubSample(cfg,EEG)
sampleSize = cfg.subsample.length;

if cfgcheck(cfg.subsample,'startpoint','random') %if no startpoint specified, pick a random one
    if sampleSize < length(EEG.data)
        startPoint = randi(length(EEG.data)-sampleSize);
    else
        startPoint = 1;
    end
else %specify a specific latency to start at
    startPoint = cfg.subsample.startpoint;
end

disp(' ')
disp(['Processing data from data point ' num2str(startPoint) ' to data point ' num2str(startPoint+sampleSize) '...'])

if startPoint+sampleSize < length(EEG.data)
    EEG.data = EEG.data(:,startPoint:(startPoint+sampleSize));
else
    warning('Not enough samples to meet sample size requirement')
    disp(' ')
    disp(['Continuing with ' num2str(length(EEG.data)-startPoint) ' samples...'])
    disp([num2str(startPoint+sampleSize-length(EEG.data)) ' samples missing...'])
    EEG.data = EEG.data(:,startPoint:end);
end
end

