function parsave(filepath,varargin)

count = 1;
for c = 1:2:length(varargin)
    varnames{count} = varargin{c};
    input = varargin{c+1};
    eval([varnames{count} '=input;']);
    count = count+1;
end
clear input varargin count c

save(filepath,varnames{:},'-v7.3')
