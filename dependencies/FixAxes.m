function FixAxes(ax,fontsize)

if ~exist('fontsize','var')
   fontsize = 14; 
end

set(ax, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'FontSize'    , fontsize        , ...
    'LineWidth'   , 1.5         , ...
    'TickLength',0.5*[0.01 0.025]);

ax.XLabel.Color = [0.15 0.15 0.15];
ax.YLabel.Color = [0.15 0.15 0.15];
set(gcf,'color','w')

% ticks = get(ax,'XTick');
% % while any(rem(ticks,1))
% %     ticks = ticks/10;
% % end
% % scalefact = get(ax,'XTick')./ticks;
% % scalefact = scalefact(1);
% ticks = unique(round(ticks,5,'significant'));
% count = 5;
% while length(ticks) > 6
%    ticks = unique(round(ticks,count,'significant'));
%    count = count - 1;
%    if count == 0 
%       break; 
%    end
% end
% set(ax,'XTick',ticks);
% 
% 
% ticks = get(ax,'YTick');
% ticks = unique(round(ticks,5,'significant'));
% count = 5;
% while length(ticks) > 6
%    ticks = unique(round(ticks,count,'significant'));
%    count = count - 1;
%    if count == 0
%       break; 
%    end
% end
% set(ax,'YTick',ticks);



