function [eulerAngles,rotCentre,DVC_new_data] = ...
                    shoemake_3D_v04_07_02_Abdo(inDir,inFile,DUunits)

%% shoemake_3D_v04_07 (alpha)
%{
%     Removes 3-D rigid body translation and rotation from a volumetric
%     dataset of displacement data (input as DAT-type format, as from
%     DaVis). In addition to data outputs for later analysis, several
%     visualisation options are available, and may be selected in Section 2
%     ("User Inputs"), including options to realign all the data to an
%     arbitrary reference frame, export the data as a DAT-type file,plot
%     orthogonal slices, bin data and plot as separate orthogonal slices
%     and plot variation in average and standard deviation in z-direction.
%     The majority of inputs for code may be entered in "User inputs"
%     section, but some modifications must be made within the code.
%
% INPUT DATA FORMAT (DAT-type)
%   Data should be in a 6 column matrix of
%        [Coordinates(X,Y,Z) Displacements(dX,dY,dZ)]
%   Coordinates should be ordered to agree with:
%       [YY XX ZZ] = meshgrid(Ymax:-1:Ymin, Xmin:Xmax, Zmin:Zmax)
%       Coordinates = [XX(:) YY(:) ZZ(:)]  
%
% OUTPUTS
%   Structures: initial     {Data, Size, DX}
%               C           {Data, Size, DX}
%               DVC         {Data, Size, DX}
%               strainAll   {Data, Size, DX}
%               bin         {Data, Size, DX}
%               binStrain   {Data, Size, DX, Comp}
%
% Code originally developed by M. Mostafavi, with additions by M. Jordan.
% Based on Ken Shoemake's Euler angle extraction.
%
% Last edit: April 2014 (M Jordan)
% (C) 2014, University of Oxford
%}

%% 1. Routine Stuff
% format compact
% clc 
% clear
% set(0,'DefaultFigureWindowStyle','docked')
close all; 
% pause off
version = 'shoemake_3D_v04.07 (alpha)';
disp(['Running ' version '.']); 
disp ' ';

%% 2. (a) User inputs: details of dataset and outputs
%Path for input file (.dat)
% inDir = ['P:\Shraddha\Displacement\450N-4500N 290k cycles_DVC_Vec_96x96x96_75ov=unknown' '\'];

% inFile = 'B00001.DAT';
cycleName = 'Corr';

% Overwrite variable "data" already in memory and repeat calculations? 
Overwrite = true;          %(true/false)

% Dataset information not recoverable from input file
res = 5;               % resolution in mm per pixel
ov = 75;                    % The percentage overlap setting in final pass

% Crop notch? (For EE8519 only) 
CropNotch = false;           %(true/false)
notchDiameter = 4;          % 4mm at -8.1d, 2mm at -6.6d

% Crop borders (removes spurious DVC values at boundaries)
CropBorders = false;         %(true/false)
borderSize = 2;             %Size (in voxel) of border effects

% Correct for rigid body motions
CorrectRBM = true;          %(true/false)

% Reference displacements to geometry? 
GeometryReference = false;   %(true/false)
forceTheta = [30 10 20];   %**NI** 3 angles along principal axes in degrees e.g. [45 -45 60] (default:0 or false)

% Bin data? 
BinData = true;            %(true/false)
numberOfBins = 5;
PlotBinnedSlice = true;        %(true/false) Requires PlotOrthoSlices

% Extract orthogonal strains?
ExtractStrain = false;      %(true/false)

% Output options (true/false)
WriteDataToFile = true;       %(true/false)
outFolder = '\Shoed';
writeRotated = true;
printHeader = false;
requestWriteAll = true;       % (true/false) as a .dat type file, as from DaVis
requestedComponents = [];   % 3D coords + {Ux,Uy,Uz = 1,2,3 respectively}

PlotOrthoSlices = true;       %(true/false) If TRUE go to 2.1

PlotAverageAndStDev = true;   %Calc & plot the x-y slice average variation     
%>>Exclude: [#slicesFromTop #slicesFrombottom]: [0 0] is entire data set
exclusionRange = [0 0];
xScale = 'pixel';                     %x-axis in 'pixel' or 'mm'

%% 2.1 Orthoslice options
% List requested plots in 4 element row vectors: 
%{
% PLOTS = {[orientation,sliceNo,map,Vectors] 'Title'}
% ORIENTATION (for map): {Slice i=sliceNo, where i=1-3 (X-Z)};
% SLICENO: valid range: 1 <= sliceNo <= lx,ly,lz (depends on orient);
Fractional values (0<=sliceNo<1) will also work for the whole volume.
% MAP: {1-3 = Ux-Uz, 4-6 = Exx-Ezz} N.B. 4-6 will require manual alteration
% of code (search 'mapData' and change variable name as appropriate);
% Vectors: true/false; 
% Placeholder 'Title' generates either Default figure title or 'string' as required;
%
% e.g. Plots = {...
%     [2 100 3 true] 'Title';...
%     [1 50 2 false] 'A';...
%     };
% generates 2 plots: 
% 1st: XZ slice (Y=100), Uz map, vectors & default title, 
% 2nd: YZ (X=50),Uy,no vectors and titled 'A' 
%}

% %Either specify particular plots ('Plots' can be any number of plots >0)
Plots = {...
    [3 0.5 2 true] 'Title';...
    [2 0.5 2 true] 'Title';...
    [1 0.5 2 true] 'Title';...
% %     [3 19 1 true] 'Title';...
%     [3 0.5 2 true] 'Title';...
%     [3 0.5 1 true] 'Title';... 
%     [2 0.5 2 true] 'Title';...
%     [1 0.5 2 true] 'Title';...
    };

% %Or for a series of plots, uncomment the following:
% first = 2; int = 1; last = 4; 
% plot_nos = first:int:last;
% Plots = repmat({[1 0 1 true] 'Title'}, size(plot_nos,2),1);
% for i= 1:size(plot_nos,2)
%     Plots{i,1}(2) = plot_nos(i);
% end

% Plot separate figure with original input data for comparison?
PlotOriginal = true;  %(true/false)

% Specify an offset origin for the plots 
%Coordinate origin options: 'Natural', 'Centre', 'Rotation' or arbitrary [x y z] in mm;
coordinateOrigin = 'Centre';
%Displacement origin options: 'Absolute', 'Relative' or arbitrary [dx dy dz] in mm;
displacementsOrigin = 'Relative';

%For map
% Smoothing function
% F = 1; 
F = fspecial('gaussian', 3, 0.5);
% Number of contours
CN = 24;
%For vector field
% Plot every nth vector on a square array or specify as [nx ny nz]
n =1;
% Arrow head size
A_S = 1;
% Vector colour
A_C = [0 0 0];
% Arrow thickness
A_T = 2;
% Arrow size
V_S = 2;
% Units for axes
units = 'voxel';   %'mm' or 'voxel'

outOpts = {F,CN,n,A_S,A_C,A_T,V_S,units};

%% End of Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3. Import data, clean and remove rigid body motion
% Data import and calculations are slow, and so will only occur if new data
% is unique or overwrite requested. END is currently just before geometry
% reference section.

if ~exist('data','var') || Overwrite || ~strcmp(datPath,[inDir '\' inFile]);
    %% 3.1 Import data (if required) and data size
    datPath = [inDir '\' inFile];
%     cd(inDir);
    disp(['Importing data for ' cycleName])
    importDat(datPath);
    
%%%
% [X, Y, Z, Ux, Uy, Uz] = V2D(data);
% plotAllDis(X,Y,Z,Ux,Uy,Uz,DUunits)
% saveas(gcf,[fileparts(datPath) '\Raw_DVC.fig']);
% saveas(gcf,[fileparts(datPath) '\Raw_DVC.tif']); close
%%%
data = data (:,1:6);
    % Extract size of data from heading line of .dat file
    %> if error thrown comment out and uncomment alternative method
    lx = str2double(colheaders{1,4}(3:end-1));    %NB curly brackets      
    ly = str2double(colheaders{1,5}(3:end-1));    %NB curly brackets
    lz = str2double(colheaders{1,6}(3:end));      %NB, runs to 'end' not 'end-1'
    
    % Alternative method:extract size of the data when no header line
%     [lx ly lz dims] = bruteSize(tdata);
            
    % Distance between interrogation points in pixels
    dx_mm = data(2,1)-data(1,1);
    dx = dx_mm/res;
    % dx = (100/(100-ov))*(data(2,1)-data(1,1))/res*10^3;
    
    % Collect current data, dataset size and voxel size together
    %> To access information stored in a struct, use e.g. C.Data or C.Size
    C = struct('Data',data,'Size',[lx ly lz],'DX', dx_mm);   %C for current
    disp(['Size of data /voxel: (lx ly lz) = (' num2str(C.Size) ').']);    
    
    % Definitions
    component = {'U_x','U_y','U_z','|U|'};
    outDir = [inDir outFolder];
    eulerAngles = [nan;nan;nan];
    notchOrigin = [0,0];
    theta_deg = [0 0 0];
    
    %% 3.2 Crop notch (only for EE8519)

    if CropNotch
        [C.Data,notchOrigin,theta_deg] = crop_notch(C.Data,notchDiameter,C.Size,dx);       
    end
    
    if forceTheta == 0
        %do nothing
    else
        theta_deg = forceTheta;
    end
    
    %% 3.3 Cleaning: Replace rows with 0 displacements to NaN
    
    tdata_x = C.Data(:,4);
    tdata_y = C.Data(:,5);
    tdata_z = C.Data(:,6);
    notn = (tdata_x==0 & tdata_y==0 & tdata_z==0);
    tdata_x(notn) = NaN;
    tdata_y(notn) = NaN;
    tdata_z(notn) = NaN;

    C.Data(:,4:6) = [tdata_x tdata_y tdata_z];
    
    %% 3.3.1 Removing border effects
    if CropBorders
        smallData = C.Data; 
        smallData = smallData(smallData(:,1)> borderSize*C.DX & smallData(:,1)<(C.Size(1)-borderSize)*C.DX,:);
        smallData = smallData(smallData(:,2)> borderSize*C.DX & smallData(:,2)<(C.Size(2)-borderSize)*C.DX,:);
        smallData = smallData(smallData(:,3)> borderSize*C.DX & smallData(:,3)<(C.Size(3)-borderSize)*C.DX,:);

        C.Data = smallData;
        C.Size = C.Size-2*borderSize;
    end
    
    initial = C;
    %% 3.3.2 Clean dataset, and correct for rigid body motions
    if  CorrectRBM
        [C.Data,eulerAngles,rotCentre] = RBMCorr('Minimum',C.Data,C.Size);
    else 
        disp 'No corrections for RBM made.'
    end
        
    %% 3.4.7 Reform to DaVis-like file (inc 0 displacement coordinates)
    if CropBorders
        C.Data(:,1:3) = smallData(:,1:3);
    else
        C.Data(:,1:3) = data(:,1:3);
    end

    DVC_new_data=C.Data;
    DVC_new_data(isnan(DVC_new_data))=0;
    
    DVC = struct('Data',DVC_new_data,'Size',C.Size,'DX', dx_mm);
    
    correctedData = C.Data;
    disp(['Size of corrected data (pre-rotation) /voxel (lx ly lz) = (' num2str(C.Size) ').']);
else
    disp('*** Using previous results until geometry reference section (S4.1) ***')
    C.Size = initial.Size;
    C.DX = initial.DX;
    disp('Dataset details:')
    disp(['   Size of corrected data /voxel: (lx ly lz) = (' num2str(C.Size) ').']);    
    disp(['   Euler angles /degrees: (alpha beta gamma) = (' num2str(radtodeg(eulerAngles),'% 1.3f ') ').'])
    disp('****')
end
disp(' ')

%% 4 Rotate coordinates and displacements to align with geometry
C.Data = correctedData;
% theta_deg=30; disp(['theta_deg overwritten (= ' num2str(theta_deg) ')'])

if GeometryReference && (mod(sum(abs(theta_deg)),360) ~= 0)
    [C.Data, C.Size, C.DX] = geoRef3D(C.Data, C.Size, C.DX, theta_deg);
end
    %{
if GeometryReference && (mod(sum(abs(theta_deg)),360) ~= 0)    
    
    [C.Data(:,4:6)] = dispNotchRef(C.Data(:,4:6), theta_deg);
% OC6 % OC8
    [C.Data, C.Size, C.DX] = imRotateShell(C.Data,C.Size,C.DX,theta_deg);
    end
%}

%% 5. Quantification
%% 5.1 Strain calculation
if ExtractStrain
    %An inelegant method of collecting Exx,Eyy and Ezz and coordinates in a
    %columnwise data format. calcOrthoStrain is embedded in plotOrthoSlice
    %(S.7) and so it is not necessary to compute the strains here if only
    %wish to plot them. If binning data, strains calculated here may be binned afterwards(S5.2),
    %while strains calculated in plotOrthoSlice will be from binned data.
    
    [dataOut,coordsOut,sizeOut,compOut] = calcOrthoStrain(4,C.Data,C.Size,'simple',C.DX);
    strainDataOut = zeros(prod(sizeOut),6);
    strainDataOut(:,1:4) = [coordsOut dataOut];
    
    for i=5:6
        [strainDataOut(:,i),coordsOut,sizeOut,compOut] = calcOrthoStrain(i,C.Data,C.Size,'simple',C.DX);
    end
    strainAll = struct('Data',strainDataOut,'Size',sizeOut,'DX',C.DX,'Component',['All']);    
end
% OC5
%% 5.2 Data binning along notch direction

clear bin binStrain
if BinData
% OC7    
    binDir = 1;
    [dataOut,sizeOut,dxOut] = xbinData(binDir,C.Data,C.Size,C.DX,numberOfBins);
    bin = struct('Data',dataOut,'Size',sizeOut,'DX', dxOut);
    disp(['Binned data available >> using ' num2str(numberOfBins) ' bins.'])
    
    if ExtractStrain
        binDir = 1;
        [dataOut,sizeOut,dxOut] = xbinData(binDir,strainAll.Data,strainAll.Size,strainAll.DX,numberOfBins);
        binStrain = struct('Data',dataOut,'Size',sizeOut,'DX', dxOut);
        disp(['Binned strain data available >> using ' num2str(numberOfBins) ' bins.'])              
    end
else disp('No binned data available')
end

%% 6. Write rotated data to file

if WriteDataToFile
    if writeRotated
        writeThis = C.Data;    
        writeSize = C.Size;
    else
        writeThis = DVC.Data;
        writeSize = DVC.Size;
    end
    
    writeThis(isnan(writeThis))=0;    
    writeDataToFile(requestedComponents,requestWriteAll,printHeader,inDir,...
        outFolder,writeThis,writeSize,cycleName,[],eulerAngles)
    clear writeThis writeSize
    % OC1
else disp 'Write data to file not requested.'
end

disp(' ')

%% 7. Plotting of 2D ortho slices
% OC2
if PlotOrthoSlices
% OC4
    if strcmp(coordinateOrigin,'Rotation')
        if exist('rotCentre','var')
            coordOrigin = rotCentre;
        else
            disp 'Rotation centre underfined. Orthoslices centred on (0,0,0)'
            coordOrigin = 'Natural';
        end
    else 
        coordOrigin = coordinateOrigin;
    end
    
    borderWidth = borderSize*CropBorders;
    
    for i=1:size(Plots,1)
        %Slice orientation, number, requested map and vector option
        mapID = Plots{i,1};  %[2 Y Map Vectors]; 
        if strcmp(Plots{i,2},'Title')
            plotTitle = {true;'(rot. corr.)'; cycleName};
        else
            plotTitle = {false;Plots{i,2}; ''};
        end
        disp(['Plot ' num2str(i) ': {[' num2str(mapID) '] ' Plots{i,2} '}...']);
        
        %Requested slice number is adjusted to account for rotation of
        %sample
        if mapID(2)<1
            mapID_temp=ceil(mapID(2)*C.Size(mapID(1)));
        else
            mapID_temp=ceil(C.Size(mapID(1))*mapID(2)/initial.Size(mapID(1)));
        end

        mapID_corr = [mapID(1) mapID_temp-borderWidth mapID(3:end)];
        disp(['** Plotting corrected slice: ' num2str(mapID_temp) ' **'])
        
        figure
%         plotOrthoSlice(mapID_corr,C.Data,C.Size,...
%             coordOrigin,displacementsOrigin,outOpts,plotTitle,C.DX);
%         caxis([-0.0084 0.0091]); disp '##caxis altered'
        if PlotOriginal
            % Plot original data for comparison
            plotTitle = {true;'(original)'; cycleName};

            %Requested slice number is adjusted to account for rotation of
            %sample
            if mapID(2)<1
                mapID_temp=ceil(mapID(2)*initial.Size(mapID(1)));
            else
                mapID_temp=mapID(2);
            end
            
            mapID_init = [mapID(1) mapID_temp-borderWidth mapID(3:end)];
            disp(['** Plotting original slice: ' num2str(mapID_temp) ' **'])            
            
%             figure
%             plotOrthoSlice(mapID_init,initial.Data,initial.Size,...
%                 coordOrigin,displacementsOrigin,outOpts,plotTitle,initial.DX);
%             caxis([0.0084 -0.0091]); disp '##caxis altered'
        end
        
        if PlotBinnedSlice && exist('bin','var')
            % Plot binned data
            plotTitle = {true;'(binned)';cycleName};
            outOptsBin = {F,CN,[1 n n],A_S,A_C,A_T,V_S,units};
            %Requested slice number is adjusted to account for rotation of
            %sample
            if mapID(2)<1
                mapID_temp=ceil(mapID(2)*bin.Size(mapID(1)));
            else
                mapID_temp=ceil(bin.Size(mapID(1))*mapID(2)/initial.Size(mapID(1)));
            end
            
            mapID_bin = [mapID(1) mapID_temp-borderWidth mapID(3:end)];
            disp(['** Plotting binned slice: ' num2str(mapID_temp) ' **'])
            figure
            plotOrthoSlice(mapID_bin,bin.Data,bin.Size,...
                coordOrigin,displacementsOrigin,outOptsBin,plotTitle,bin.DX);
%             caxis([-0.0084 0.0091]); disp '##caxis altered'
            % OC9
        end
    end
else disp 'Ortho slices not requested.'
end

disp(' ')
%% 8. Plotting Average and StDev (slice-wise) in Z-direction

if PlotAverageAndStDev 
    z_range = [exclusionRange(1)+1 C.Size(3)-exclusionRange(2)];
    z_length = z_range(2)-z_range(1)+1;

    % RMS and Ave have column headings: 
    %
    % ||Slice_number|Z_coord|Blank!|<Ux>|<Uy>|<Uz>|<U>|{SD(U)/<U>|}|

    [RMS, ave] = stdevAndAverage(C.Data(:,4:6),C.Size,z_range,C.DX);   
    %OC3
    figure;  
    hold on
    switch xScale
        case 'pixel'
            xlabel('Distance along Z /pixel');
            xlim(z_range);
            xCol = 1;
        case 'mm'
            xlabel('Distance along Z /mm');
            xlim(z_range*dx_mm);
            xCol = 2;
        otherwise
            disp('ERROR: Invalid X-axis unit; valid units are {pixel/mm}')
            ERROR
    end
    %     plot(Ave(:,xScale),Ave(:,4),'rx-',Ave(:,xScale),Ave(:,5),'gx-',Ave(:,xScale),Ave(:,6),'bx-',Ave(:,1),Ave(:,7),'k.-',[0,max(RMS(:,1))],[0,0],'k-');
    errorbar(ave(:,xCol),ave(:,4),RMS(:,4)./2,'rx-');
    errorbar(ave(:,xCol),ave(:,5),RMS(:,5)./2,'gx-');
    errorbar(ave(:,xCol),ave(:,6),RMS(:,6)./2,'bx-');
    errorbar(ave(:,xCol),ave(:,7),RMS(:,7)./2,'k.-');
    plot([0,max(RMS(:,1))],[0,0],'k-');          %plot x-axis line at y=0

    title(['Average(U_i) over X-Y slices for ' cycleName]);
    ylabel('X-Y Slice Average /mm');
    legend('U_x','U_y','U_z','|U|','Location','EastOutside');
    hold off
    
else disp 'Average and StDev plots not requested.'
end

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  savefig(FigHandle, fullfile(outDir, [num2str(iFig), '.fig']));  
  saveas(FigHandle, fullfile(outDir, [num2str(iFig), '.tif'])); close
end

close all; disp(['*** ' version ' complete ***'])
% [X, Y, Z, Ux, Uy, Uz] = V2D(DVC_new_data);
% plotAllDis(X,Y,Z,Ux,Uy,Uz,DUunits)
% saveas(gcf,[inDir '\Corr_DVC.fig']);
% saveas(gcf,[inDir '\Corr_DVC.tif']); close
% end