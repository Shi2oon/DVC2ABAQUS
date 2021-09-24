function [files,Operation,unit,stif] = Westergaard_3D(StressIntensityFactor,Mode,stPs,file2)
% close all; clear; clc
% check mode   
if ~exist('Mode','var');        Mode = 'I';         end
if isempty(Mode);               Mode = 'I';         end

if exist('file2','var')
    newfile = fullfile(file2, [num2str(StressIntensityFactor) '_' Mode ' WS']);
else
    newfile = [ pwd '\3D_' num2str(StressIntensityFactor) '_' Mode ' WS']; 
end
mkdir(newfile);


minGrid = -4E-3; % m
if exist('stPs','var')
    gridStep = 8e-3/stPs;
else
    gridStep = 0.2E-3; % m
end

maxGrid = 4E-3; % m

vec = minGrid : gridStep : maxGrid; % m
[x,y,z] = meshgrid(vec,vec,vec); % m
% StressIntensityFactor = 30; % [MPa m^0.5]
fprintf('preparing synthetic Westergaard Solution Data .. ');
K = StressIntensityFactor * 1E6; % Stress intensity factor [Pa m^0.5]
E = 210E9; % Young's Modulus [Pa]
nu = 0.3; % poisson ratio
mu = E/(2.*(1+nu)); % Shear Modulus [Pa]

kappa = 3 - (4 .* nu); % [/]

% theta r and are polar coordinates centered at the crack tip in a plane 
% orthogonal to the crack front.
[theta,r] = cart2pol(x,y,z);
switch Mode
    case 'I'
        ux = (K./(2.*mu)).*sqrt(r/(2*pi)).*cos(theta/2).*(kappa - 1 + ...
              2.*(sin(theta/2)).^2); % Anderson p99 A2.44a
        uy = (K./(2.*mu)).*sqrt(r/(2*pi)).*sin(theta/2).*(kappa + 1 - ...
              2.*(cos(theta/2)).^2); % Anderson p99 A2.44b
        uz = zeros(size(uy));
    case 'II'
        ux = (K./(2.*mu)).*sqrt(r/(2*pi)).*sin(theta/2).*(kappa + 1 + ...
              2.*(cos(theta/2)).^2); % Anderson p99 A2.44a
        uy = (K./(2.*mu)).*sqrt(r/(2*pi)).*cos(theta/2).*(kappa - 1 - ...
              2.*(sin(theta/2)).^2); % Anderson p99 A2.44b
        uz = zeros(size(uy));
    case 'III'
        uz = (2*K./mu).*sqrt(r/(2*pi)).*sin(theta/2); % Anderson p99 A2.44b
        ux = zeros(size(uz)); % Anderson p99 A2.44a
        uy = zeros(size(uz)); % Anderson p99 A2.44a
    case 'fun'
        ux = ((K./(2.*mu)).*sqrt(r/(2*pi)).*cos(theta/2).*(kappa - 1 + 2.*(sin(theta/2)).^2)...
             +(K./(2.*mu)).*sqrt(r/(2*pi)).*sin(theta/2).*(kappa + 1 + 2.*(cos(theta/2)).^2))./3; 
        uy = ((K./(2.*mu)).*sqrt(r/(2*pi)).*sin(theta/2).*(kappa + 1 -  2.*(cos(theta/2)).^2)...
             +(K./(2.*mu)).*sqrt(r/(2*pi)).*cos(theta/2).*(kappa - 1 - 2.*(sin(theta/2)).^2))./3;  
        uz = (2*K./mu).*sqrt(r/(2*pi)).*sin(theta/2)./3;
        
end
plotAllDis(x,y,z,ux,uy,uz,'m')
saveas(gcf, [newfile '\' Mode '_Disp_fields.tiff']);   
saveas(gcf, [newfile '\' Mode '_Disp_fields.fig']); close  

Plot3D(sqrt(ux.^2+uy.^2+uz.^2),x,y,z,'m','U_{mag}')
saveas(gcf, [newfile '\' Mode '_Disp_Mag.tiff']);   
saveas(gcf, [newfile '\' Mode '_Disp_Mag.fig']); close  
     
    %% save displacement fields
    Operation{1,1} = 'DVC';               unit{1,1}  = 'm';
    alldata = [x(:) y(:) z(:) ux(:) uy(:) uz(:)]; % m
    if exist('stPs','var')
        files{1,1} = [newfile '\' num2str(StressIntensityFactor) 'MPa_m_DISP_' num2str(stPs) '.mat'];
    else
        files{1,1} = [newfile '\' num2str(StressIntensityFactor) 'MPa_m_DISP.mat'];
    end
    save(files{1,1},'alldata')
    
    Operation{1,2} = 'DVC';               unit{1,2}  = 'mm';
    alldata = [x(:) y(:) z(:) ux(:) uy(:) uz(:)].*1e3; % mm
    if exist('stPs','var')
        files{1,2} = [newfile '\' num2str(StressIntensityFactor) 'MPa_mm_DISP_' num2str(stPs) '.mat'];
    else
        files{1,2} = [newfile '\' num2str(StressIntensityFactor) 'MPa_mm_DISP.mat'];
    end
    save(files{1,2},'alldata')
    
    Operation{1,3} = 'DVC';               unit{1,3}  = 'um';
    alldata = [x(:) y(:) z(:) ux(:) uy(:) uz(:)].*1e6; % um
    if exist('stPs','var')
        files{1,3} = [newfile '\' num2str(StressIntensityFactor) 'MPa_um_DISP_' num2str(stPs) '.mat'];
    else
        files{1,3} = [newfile '\' num2str(StressIntensityFactor) 'MPa_um_DISP.mat'];
    end
    save(files{1,3},'alldata')

    %% save strain fields
    [exx,dux_dy,dux_dz] = gradient(ux,gridStep);
    [duy_dx,eyy,duy_dz] = gradient(uy,gridStep);
    [duz_dx,duz_dy,ezz] = gradient(uz,gridStep);
    exy = 0.5*(dux_dy + duy_dx);
    exz = 0.5*(dux_dz + duz_dx);
    eyz = 0.5*(duz_dy + duy_dz);
    o.Exy = exy;    o.Exz = exz;    o.Eyz = eyz;    o.Ezz = ezz;
    plotAllDis(x,y,z,exx,eyy,o,'abs');
    saveas(gcf, [newfile '\' Mode '_Strain_fields.tiff']);   
    saveas(gcf, [newfile '\' Mode '_Strain_fields.fig']); close  
    
    Operation{2,1} = 'Str';               unit{2,1}  = 'm';
    alldata = [x(:) y(:) z(:) exx(:) eyy(:) ezz(:) exy(:) exz(:) eyz(:)]; % m
    if exist('stPs','var')
        files{2,1} = [newfile '\' num2str(StressIntensityFactor) 'MPa_m_Strain_' num2str(stPs) '.mat'];
    else
        files{2,1} = [newfile '\' num2str(StressIntensityFactor) 'MPa_m_Strain.mat'];
    end
    save(files{2,1},'alldata')
    
    Operation{2,2} = 'Str';               unit{2,2}  = 'mm';
    alldata = [x(:).*1e3 y(:).*1e3 z(:).*1e3 exx(:) eyy(:) ezz(:) exy(:) exz(:) eyz(:)]; % mm
    if exist('stPs','var')
        files{2,2} = [newfile '\' num2str(StressIntensityFactor) 'MPa_mm_Strain_' num2str(stPs) '.mat'];
    else
        files{2,2} = [newfile '\' num2str(StressIntensityFactor) 'MPa_mm_Strain.mat'];
    end
    save(files{2,2},'alldata')
    
    Operation{2,3} = 'Str';               unit{2,3}  = 'um';
    alldata = [x(:).*1e6 y(:).*1e6 z(:).*1e6 exx(:) eyy(:) ezz(:) exy(:) exz(:) eyz(:)]; % um
    if exist('stPs','var')
        files{2,3} = [newfile '\' num2str(StressIntensityFactor) 'MPa_um_Strain_' num2str(stPs) '.mat'];
    else
        files{2,3} = [newfile '\' num2str(StressIntensityFactor) 'MPa_um_Strain.mat'];
    end
    save(files{2,3},'alldata')
    
    stif{1}='E';
    stif{2}='E';
    stif{3}='A';
    
    %%
e33 = zeros(size(exx(:,:,1)));
Maps.S11_F = e33;       Maps.S22_F = e33;       Maps.S12_F = e33;
Maps.S21_F = e33;       Maps.S31_F = e33;       Maps.S13_F = e33;
Maps.S32_F = e33;     	Maps.S23_F = e33;       Maps.S33_F = e33;

    Data.X = x(:,:,1) *1e6;       	Data.Y = y(:,:,1) *1e6; %um 
    Maps.E11_F = exx(:,:,1);	Maps.E12_F = exy(:,:,1); 	Maps.E13_F = exz(:,:,1);
    Maps.E21_F = exy(:,:,1); 	Maps.E22_F = eyy(:,:,1);  	Maps.E23_F = eyz(:,:,1);
    Maps.E31_F = exz(:,:,1); 	Maps.E32_F = eyz(:,:,1);  	Maps.E33_F = ezz(:,:,1);  
    Maps.W11_F1= e33(:,:,1); 	Maps.W12_F1= e33(:,:,1);   	Maps.W13_F1= e33(:,:,1);
    Maps.W21_F1= e33(:,:,1);  	Maps.W22_F1= e33(:,:,1);  	Maps.W23_F1= e33(:,:,1);
    Maps.W31_F1= e33(:,:,1);  	Maps.W32_F1= e33(:,:,1); 	Maps.W33_F1= e33(:,:,1);
    
    Maps.PH_2 = e33(:,:,1);   	Maps.MAE_2=e33(:,:,1);    	Map_RefID = e33(:,:,1); 
    Maps.E  = E/1e9; %[GPa]
    iPut.stiffnessvalues = [230 135 135 0 0 0; 135, 230, 135, 0 0 0; 135 135 230 0 0 0;...
                      0 0 0 117 0 0; 0 0 0 0 117 0 ; 0 0 0 0 0 117 ];         
    Maps.nu = nu;
    Data_InputMap.X_axis = unique(Data.X); Data_InputMap.Y_axis = unique(Data.Y);
    GrainData.RefPoint.x=0;GrainData.RefPoint.y=0;
    
    MicroscopeData.NROWS = size(ezz,1); MicroscopeData.NCOLS = size(ezz,1); 
    files{3,3} = [newfile '\' num2str(StressIntensityFactor) 'MPa_um_XEBSD.mat'];
    Operation{3,3} = 'xED';         unit{3,3}  = 'um';
    save(files{3,3},'Map_RefID','Maps','Data','iPut','GrainData','Data_InputMap',...
        'MicroscopeData')
fprintf ('DONE\n\n');