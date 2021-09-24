function plotAllDis(dvcdatapth,units,MainDispl)
dvcdata = importdata(dvcdatapth);
dvcdata = dvcdata.data;
X = reshape(dvcdata(:,1),length(unique(dvcdata(:,1))),length(unique(dvcdata(:,2))),length(unique(dvcdata(:,3))));
Y = reshape(dvcdata(:,2),length(unique(dvcdata(:,1))),length(unique(dvcdata(:,2))),length(unique(dvcdata(:,3))));
Z = reshape(dvcdata(:,3),length(unique(dvcdata(:,1))),length(unique(dvcdata(:,2))),length(unique(dvcdata(:,3))));
Ux = reshape(dvcdata(:,4),length(unique(dvcdata(:,1))),length(unique(dvcdata(:,2))),length(unique(dvcdata(:,3))));
Uy = reshape(dvcdata(:,5),length(unique(dvcdata(:,1))),length(unique(dvcdata(:,2))),length(unique(dvcdata(:,3))));
Uz = reshape(dvcdata(:,6),length(unique(dvcdata(:,1))),length(unique(dvcdata(:,2))),length(unique(dvcdata(:,3))));

if exist('MainDispl','var')==0;     MainDispl = 'Uy';   end
if exist('units','var')==0;         units = 'mm';       end

close all;    
h1 = subplot(1,3,1); Plot3D(Ux,X,Y,Z,units,'U_x');      colorbar off; 
h2 = subplot(1,3,2); Plot3D(Uy,X,Y,Z,units,'U_y');      colorbar off; 
h3 = subplot(1,3,3); Plot3D(Uz,X,Y,Z,units,'U_z');      colorbar off;  
cbax  = axes('visible', 'off');             
eval(sprintf('caxis(cbax,[min(%s(:)) max(%s(:))]);',MainDispl,MainDispl));
h = colorbar(cbax, 'location', 'southoutside','position', [0.3513 0.0993 0.3 0.03] );
h.Label.String = ['\DeltaU [' units ']']; 
set([h1 h2 h3],"clim",caxis);           set(gcf,'position',[1 41 1920 963]);  
end

function Plot3D(M,X,Y,Z,units,tiz)
if isempty(X)
    [X,Y,Z] = meshgrid(1:size(M,2),1:size(M,1),1:size(M,3)); 
    units = 'Pixel';
    hs = slice(X,Y,Z,M,unique(X),unique(Y),unique(Z)) ;
else
    hs = slice(Y,X,Z,M,unique(X),unique(Y),unique(Z)) ;
end
shading interp;                             set(hs,'FaceAlpha',1);
set(gcf,'position',[30 50 1300 950]);       
y = unique(Y);      xlim([min(y) max(y)]);
x = unique(X);      ylim([min(x) max(x)]);
xlabel(['X [' units ']']);                
ylabel(['Y [' units ']']); 
zlabel(['Z [' units ']']);
axis image; 
view([45 45 15])
c = colorbar;       colormap jet;       c.Label.String = units;
title(tiz)
end