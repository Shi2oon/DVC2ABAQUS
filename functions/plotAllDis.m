function plotAllDis(X,Y,Z,Ux,Uy,Uz,units)
figure; iV=1;
if strcmpi(units, 'abs')
    uHead = '\epsilon';
    iV=2;
    h4 = subplot(iV,3,4); Plot3D(Uz.Exy,X,Y,Z,units,[uHead '_{xy}']);      
    c  = colorbar;	cU(4,:) = c.Limits;     colorbar off;   axis off; axis tight
    h5 = subplot(iV,3,5); Plot3D(Uz.Exz,X,Y,Z,units,[uHead '_{xz}']);
    c  = colorbar;	cU(5,:) = c.Limits;     colorbar off;   axis off; axis tight
    h6 = subplot(iV,3,6); Plot3D(Uz.Eyz,X,Y,Z,units,[uHead '_{yz}']);
    c  = colorbar;	cU(6,:) = c.Limits;     colorbar off;   axis off; axis tight
    Uz = Uz.Ezz;
else
    uHead = 'U';
end
    
h1 = subplot(iV,3,1); Plot3D(Ux,X,Y,Z,units,[uHead '_x']);      
c  =colorbar;	cU(1,:) = c.Limits;     colorbar off;  axis tight
h2 = subplot(iV,3,2); Plot3D(Uy,X,Y,Z,units,[uHead '_y']);
c  =colorbar;	cU(2,:) = c.Limits;     colorbar off;  axis tight
h3 = subplot(iV,3,3); Plot3D(Uz,X,Y,Z,units,[uHead '_z']);
c  =colorbar;	cU(3,:) = c.Limits;     colorbar off; axis tight
if strcmpi(units, 'abs')
    subplot(iV,3,1); axis off; subplot(iV,3,2); axis off; subplot(iV,3,3); axis off;
end
cbax  = axes('visible', 'off');   
cU(abs(cU)==1)=0;
caxis(cbax,[min(cU(:)) max(cU(:))]);
if strcmpi(units, 'abs')
    h = colorbar(cbax, 'location', 'westoutside','position', [0.1122 0.2728 0.0143 0.4720] );
else
    h = colorbar(cbax, 'location', 'southoutside','position', [0.3513 0.0993 0.3 0.03] );
end
h.Label.String = [uHead ' [' units ']']; 
set([h1 h2 h3],"clim",caxis);         
if strcmpi(units, 'abs');   set([h1 h2 h3 h4 h5 h6],"clim",caxis);end
set(gcf,'position',[1 41 1920 1003]);  
end