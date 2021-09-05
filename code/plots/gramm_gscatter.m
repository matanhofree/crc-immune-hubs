function outfig = gramm_gscatter(x,y,grpType,clrMap,sym,siz,doleg,xnam,ynam)
%%
    gfig = gramm('x',x,'y',y,'color',grpType);
    gfig.geom_point('dodge',20);
    gfig.set_names('x','tSNE 1','y','tSNE 2','color','Type');
    if exist('clrMap','var') && ~isempty(clrMap) 
        gfig.set_color_options('map',clrMap);
    end
    
    if exist('siz','var') && ~isempty(siz) 
        gfig.set_point_options('base_size',siz);
    end
    
    % gfig.set_layout_options('position',[0.05 0.05 0.75 0.9],'legend_position',[ 0.89 0.1 0.08 0.9 ]);    
    
    gfig.draw()
    
    
    
    

end