function [s, cfg] = ft_statfun_signrank(cfg, dat, design)

tmp = unique(design);
for c = 1:length(tmp)
indices{c} = find(design == tmp(c));
orgdat(:,:,c) = dat(:,indices{c});
end

p = ones(size(dat,1),1);
for c = 1:size(dat,1)
    p(c) = friedman(squeeze(orgdat(c,:,:)),1,'off');
%     p(c) = signrank(dat(c,indices{1}),dat(c,indices{2}));
%     if median(dat(c,indices{1})) > median(dat(c,indices{2}))
%         p(c) = 1-p(c);
%     else
%         p(c) = -1+p(c);
%     end
    
end
s.stat = p;
% 
% switch cfg.tail
%     case 0
%         s.critval = [-0.95 0.95];
%     case 1
        s.critval = 0.95;
%     case -1
%         s.critval = -0.95;
% end

s.df = size(dat,2)-2;

end