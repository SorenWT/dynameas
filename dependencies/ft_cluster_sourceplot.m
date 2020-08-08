function figaxes = ft_cluster_sourceplot(datain,sourcemodel,atlas,roimask,varargin)

for c = 1:length(sourcemodel.pos)
    if isfield(atlas,'parcellation')
        plotdata(c) = datain(atlas.parcellation(c));
        plotmask(c) = roimask(atlas.parcellation(c));
    elseif isfield(atlas,'parcels')
        plotdata(c) = datain(atlas.parcels(c));
        plotmask(c) = roimask(atlas.parcels(c));
    end
end
plotmask = double(plotmask);

argsin = varargin;

argsin = setdefault(argsin,'method','brains4');

if EasyParse(argsin,'method','wholebrain')
    bnd.pnt = sourcemodel.pos;
    bnd.tri = sourcemodel.tri;
    
    ft_plot_mesh(bnd,'facealpha',vert(1-plotmask));
    ft_plot_mesh(bnd, 'vertexcolor', vert(plotdata), 'facealpha', vert(plotmask), 'maskstyle', 'opacity');
    
    lighting gouraud
    camlight
elseif EasyParse(argsin,'method','brains4')
    figaxes = brains4(datain,sourcemodel,atlas,'mask',roimask);
end

    



