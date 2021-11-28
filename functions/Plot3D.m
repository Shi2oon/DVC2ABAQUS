%% one plot
function Plot3D(M,X,Y,Z,units,tiz)
if ~exist('X','var')
    [X,Y,Z] = meshgrid(1:size(M,1),1:size(M,2),1:size(M,3)); 
    units = 'pixel';
    tiz = '';
end
try
    hs = slice(X,Y,Z,M,unique(X),unique(Y),unique(Z)) ;
    set(hs,'FaceAlpha',1);
catch
    try
        hs = slice(Y,X,Z,M,unique(Y),unique(X),unique(Z));
        set(hs,'FaceAlpha',1);
    catch
        scatter3(X(:),Y(:),Z(:),[],M(:),'filled');
    end
end
shading interp;                             
set(gcf,'position',[30 50 1300 950]);       axis image; 
y = unique(Y);      xlim([min(y) max(y)]);
x = unique(X);      ylim([min(x) max(x)]);
ylabel(['X [' units ']']);                
xlabel(['Y [' units ']']); 
zlabel(['Z [' units ']']);
c = colorbar;       colormap jet;       c.Label.String = units;
title(tiz)
end