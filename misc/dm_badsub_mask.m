function outputs = dm_badsub_mask(outputs,badsubs)
% dm_badsub_mask adds bad subject information to a Dynameas structure
%
% Inputs:
%    outputs: a Dynameas outputs structure. This can be a cell array, in
%       which case badsubs should be a cell array of equal length
%    badsubs: a list of bad subjects. This can either be a vector of
%       indices or a cell array of filenames
%
% Outputs:
%    outputs: the Dynameas output structure with the bad subject mask and
%       information added

if iscell(outputs)
    for i = 1:length(outputs)
        outputs{i} = dm_badsub_mask(outputs{i},badsubs{i});
    end
else
    
    if iscell(badsubs)
        badsub_ind = contains(outputs.sub,badsubs);
    else
        if length(unique(badsubs))==2
            badsub_ind = badsubs;
        else
            badsub_ind = unfind(badsubs,size(outputs.data,1));
        end
    end
    
    badsub_ind = vert(badsub_ind);
    
    outputs.badsub_mask = ~repmat(badsub_ind,1,size(outputs.data,2),size(outputs.data,3));
        
    outputs.data_orig = outputs.data;
    
    outputs.data_masked = outputs.data.*nanmask(outputs.badsub_mask);
    
    outputs.data = outputs.data_masked;
        
    outputs.badsubs = badsub_ind;
    
end