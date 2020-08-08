function [s, cfg] = ft_statfun_ranksum(cfg, dat, design)

tmp = unique(design);
indices{1} = find(design == tmp(1));
indices{2} = find(design == tmp(2));

switch cfg.tail
    case 0
        tail = 'both';
        s.critval = [-1+cfg.alpha 1-cfg.alpha];
    case 1
        tail = 'right';
        s.critval = 1-cfg.alpha;
    case -1
        tail = 'left';
        s.critval = -1+cfg.alpha;
end

p = ones(size(dat,1),1);
for c = 1:size(dat,1)
    p(c) = ranksum(dat(c,indices{1}),dat(c,indices{2}),'tail',tail);
    if median(dat(c,indices{1})) > median(dat(c,indices{2}))
        p(c) = 1-p(c);
    else
        p(c) = -1+p(c);
    end
    
end
s.stat = p;

s.df = size(dat,2)-2;

end