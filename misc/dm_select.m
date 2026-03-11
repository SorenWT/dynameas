function outputs = dm_select(outputs,field,selection)

if iscell(outputs)
    for i = 1:length(outputs)
        outputs{i} = dm_select(outputs{i},field,selection);
    end
else
    
    
    if isnumeric(outputs.(field))
        ind = outputs.(field)==selection;
    elseif iscell(outputs.(field)) && iscell(selection)
        ind = contains(outputs.(field),selection);
    else
        if length(unique(selection))==2
            ind = find(selection);
        end
        
        ind = selection;
    end
    
    dimtok = tokenize(outputs.dimord,'_');
    
    whichdim = find(strcmpi(field,dimtok));
    
    indmat = repmat({':'},1,length(dimtok));
    
    indmat{whichdim} = ind;
    
    outputs.data = outputs.data(indmat{:});
    outputs.(field) = outputs.(field)(ind);
    
    if isfield(outputs,'badsubs')
       outputs.data_orig = outputs.data_orig(indmat{:});
       outputs.data_masked = outputs.data_masked(indmat{:});
       outputs.badsub_mask = outputs.badsub_mask(indmat{:});

       if strcmpi(field,'sub')
           outputs.badsubs = outputs.badsubs(ind); 
       end
    end
    
end

