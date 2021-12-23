 function dispOrigin = plotOrthoSlice(mapID,data3D_col,dataSize,coordinateOrigin,displacementsOrigin,outOpts,Title,dx_mm)
% plotOrthoSlice
%     plotOrthoSlice(mapID,data3D_col,dataSize,coordinateOrigin,displacementsOrigin,outOpts,Title,dx_mm)
%       Plots a 2D slice of 3D gridded data as background map, with
%       overlaid arrows optional. Slice must be perpendicular to reference
%       axes.
%       
% mapID:    
%         a 4 element vector as mapID = [ORIENT SLICENO MAP VECTORS]
%         ORIENT (for map): {Slice i=sliceNo, where i=1-3 (X-Z)};
%         SLICENO: valid range: 1<=sliceNo<=lx,ly,lz (depends on orient);
%         MAP: {1-3 = Ux-Uz, 4-6 = Exx-Ezz} N.B. 4-6 will require manual 
%             alteration of code (search 'mapData' and change variable name 
%             as appropriate);
%         VECTORS: true/false; 
% data3D_col:
%         a 6 column matrix of data as [X Y Z DX DY DZ].
% dataSize: 
%         a 3 element vector of the dimensions of the 3D cube of data
%           as [lx ly lz].
% coordinateOrigin:
%         The origin coordinates of the background map. Specify as:
%         'Natural' - using input data native [0 0 0] 
%         'Centre' - the coordinate of the centroid (average)
%         [arbritrary_3D_coordinate] - a 3 element vector of [x0 y0 z0]
% displacementsOrigin:
%         The rigid body offset in the displacements. Specify as
%         'Absolute'
%         'Relative' - Relative to the average displacement
%         'RelToOri' - relative to the displacment at the coordinateOrigin
%         [arbritrary_3D_vector] - a 3 element vector of [dx0 dy0 dz0]
% outOpts:
%         Cell of options related to the appearance of the figures:
%             outOpts = {F,CN,n,A_S,A_C,A_T,V_S,units,PlotBorders};
%         where (defaults in <>):
%             FOR MAP:
%             F: Smoothing function: <fspecial('gaussian', 3, 0.5)>
%             CN: Number of contours in map <30>  
%             FOR VECTOR FIELD:
%             n: Plot every nth vector on square array(<10>)
%                 or specify as [nx ny nz] for non-square array
%             A_S: Arrow head size (scaling factor) <1>
%             A_C: Arrow colour as 3 element vector <[1 0 1]>
%             A_T: Arrow thickness <1>
%             V_S: Arrow size (scaling factor) <2>
%             GENERAL:
%             unit: Units for axes as 'mm' or 'voxel' <'mm'>
%             PlotBorders: Plots bounding box around all coordinates <true>      
%         N.B. enter as [] for defaults
% Title:
%   Plots have a default title form which defines the slice plane, and
%   user user defined data name. Also defines the colorbar title. Default
%   may be overwritten.
%       e.g. of default title: "X-Y slice (13 of 53) for Loaded Data"
%            of colourbar title: "U_y (corrected data) /mm"
%     Title = {Use_Default,Data_name(Title),Colourbar_Fragment}
%   where
%       Use_Default: <true>/false 
%       Data_name: A short string to describe data ("Loaded Data" above) 
%       Data_name: (if Use_Default = 0) Arbritrary string as figure title 
%       Colourbar_Fragment: Short string that appears in colourbar 
%             ("(corrected data)" above)
% dx_mm:
%   The spacing (resolution) of the data. Can be either 1 value for a cubic
%   array data, or a 3 element vector for a non-cubic array. <1>
%
% % Code outlined by M. Mostafavi.
% % Adapted to function and heavily modified by MSL Jordan.
% % (c) University of Oxford, 2014. Last modified Sept 2014
 
%% 1.Definitions and input check

if size(data3D_col,2)~=6
    disp(' ')
    disp('ERROR: data3D_col is invalid size. Valid matrix contains 6 columns: [coordData(:,1:3) vectorData(:,4:6)]'); ERROR
end

if isempty(outOpts)
    %No inputs,so use defaults
    F = fspecial('gaussian', 3, 0.5);     %Smoothing function
%     F=1;      %no smoothing
    CN = 30;                % Number of contours
    %For vector field
    n = 10;                 % Plot every nth vector on a square array
    A_S = 1;                % Arrow head size
    A_C = [1,0,1];          % Vector colour
    A_T = 1;                % Arrow thickness
    V_S = 2;                % Arrow size
    unit = 'mm';
    PlotBorders = true;
%     Smooth = true;    %redundant
else
    if size(outOpts,2)<10
        %DO NOT DELETE - change your inputs!
%        disp '##WARNING: plotOrthoSlice modified 30 Oct 2014'
%        disp '~~ outOpts{9} = (true/false) - see Matthew'
%        disp '~~~ Adds a bounding box to the orthoslice'
%        disp '~~ outOpts{10} = (true/false) - see Matthew'              %redundant
%        disp '~~~ Smoothes data before plotting. Previously applied automatically'
%        disp '(message will be removed after 21 Nov 2014)'
       outOpts = [outOpts{1:8} {1}];% {1}]; %redundant
    end
    F = outOpts{1};
    CN = outOpts{2};
    n = outOpts{3};
    A_S = outOpts{4};
    A_C = outOpts{5};
    A_T = outOpts{6};
    V_S = outOpts{7};
    unit = outOpts{8};
    PlotBorders = outOpts{9};
%     Smooth = outOpts{10}; %redundant
end

switch lower(unit)
    case 'voxel'
        uMult = repmat(dx_mm.^-1,1,6/size(dx_mm,2));
    otherwise
        uMult = ones(1,6);
end  

data3D_col = bsxfun(@times,data3D_col,uMult);
component = {['U_x /' unit];['U_y /' unit];['U_z /' unit];'Exx ';'Eyy ';'Ezz '}; %N.B. all components should be > 3 characters

%Define 3-D vector spacing if necessary
switch size(n,2)
    case 0
        disp(' ')
        disp('WARNING: Vector spacing (n) not entered; default of n = 10')
        N=[10 10 10];
    case 1
        N=[n n n];
    case 3
        N=n;
    otherwise
        disp('ERROR: Invalid number elements in n. Valid lengths are 1 or 3')
        ERROR
end

% Define titles and data for required slice
switch mapID(1)
    case 1
        x_temp = 2;
        y_temp = 3;
        vecElements = {mapID(2),1:N(2):dataSize(2),1:N(3):dataSize(3)};
        figX = ['y /' unit];
        figY = ['z /' unit];
        Tit1 = 'Y-Z plane (X = ';
        Tit2 = dataSize(1);
    case 2
        x_temp = 1;
        y_temp = 3;
        vecElements = {1:N(1):dataSize(1),mapID(2),1:N(3):dataSize(3)};
        figX = ['x /' unit];
        figY = ['z /' unit];
        Tit1 = 'X-Z plane (Y = ';
        Tit2 = dataSize(2);
    case 3
        x_temp = 1;
        y_temp = 2; 
        vecElements = {1:N(1):dataSize(1),1:N(2):dataSize(2),mapID(2)};
        figX = ['x /' unit];
        figY = ['y /' unit];
        Tit1 = 'X-Y plane (Z = ';
        Tit2 = dataSize(3);
    otherwise
        disp('ERROR: Slice orientation value (mapID(1)) incorrect. Valid values are 1-3.')
        ERROR
end

if Title{1}
%     disp '##WARNING: plotOrthoSlice.m "Title" definition changed 25/9/14'
%     disp '~~ Need to define as:'
%     disp '~~ Title = {Use_Default,Data_name(Title),Colourbar_Fragment}'
    figTitle = [Tit1 num2str(mapID(2)) ' of ' num2str(Tit2) ') for ' Title{2}];
else
    figTitle = Title{2};
end

%Define plotting origins
switch logical(true)
    case strcmpi(coordinateOrigin,'natural')
        coordOrigin = [0 0 0];
    case strcmpi(coordinateOrigin,'centre')
        coordOrigin = nanmean(data3D_col(:,1:3)).*uMult(1:3);
    case isnumeric(coordinateOrigin) && isrow(coordinateOrigin) && size(coordinateOrigin,2) == 3
        coordOrigin = coordinateOrigin.*uMult(1:3);
    otherwise
        error('Invalid coordinateOrigin.')
end

switch logical(true)
    case strcmpi(displacementsOrigin,'absolute')
        dispOrigin = zeros(1,6);
    case strcmpi(displacementsOrigin,'relative')
        dispOrigin = [nanmean(data3D_col(:,4:6)) 0 0 0].*uMult;    %NB Applies to displacement origin of map data [Ux Uy Uz Exx Eyy Ezz]
    case strcmpi(displacementsOrigin,'reltoori')
        %using above coordOrigin (in voxel)#
        if ~strcmpi(unit,'voxel')
            cOri = coordOrigin./dx_mm;
        else
            cOri = coordOrigin;
        end
        %Get related displacement
        cOri = round(cOri)+[1 1 1];
        oriInd = sub2ind(dataSize,cOri(1),dataSize(2)-cOri(2),cOri(3));
%         oriInd = sub2ind(dataSize,cOri(1),cOri(2),cOri(3));
        oriDisp = data3D_col(oriInd,4:6);         %data_3D_col is in output unit
        if isnan(sum(oriDisp))
            oriDisp = [0 0 0];
        end
        dispOrigin = [oriDisp 0 0 0];
    case isnumeric(displacementsOrigin) && isrow(displacementsOrigin) && size(displacementsOrigin,2) == 3
        dispOrigin = [displacementsOrigin 0 0 0].*uMult;
    otherwise
        error('Invalid displacmentsOrigin.')
end

c0 = [coordOrigin(x_temp) coordOrigin(y_temp)];
d0 = [dispOrigin(x_temp) dispOrigin(y_temp)]; %only applies to vectors

%% 2. Extract map and vector coordinates and values for given slice

% Values for map and coordinates:n = 1-3 are displacements (vector components), 4-6 are
% strains (previously calculated)
switch mapID(3)  
    case {1 2 3}
        % Values for map
        mapThis = reshape(data3D_col(:,mapID(3)+3),dataSize);

        switch mapID(1)
            case 1
                mapElements = {mapID(2),1:dataSize(2),1:dataSize(3)};
            case 2
                mapElements = {1:dataSize(1),mapID(2),1:dataSize(3)};
            case 3
                mapElements = {1:dataSize(1),1:dataSize(2),mapID(2)};
        end
        % Coordinates for map
        xx_temp=reshape(data3D_col(:,x_temp),dataSize);
        yy_temp=reshape(data3D_col(:,y_temp),dataSize);
        ssX=(squeeze(xx_temp(mapElements{:}))-c0(1));
        ssY=(squeeze(yy_temp(mapElements{:}))-c0(2));
        
    case {4 5 6}
        disp('** Strain map requested; calculating using simple method **')
        % Values for map
        [sData,sCo,sSize,sComp] = calcOrthoStrain(mapID(3),data3D_col,dataSize,'Simple',dx_mm);
        strain = struct('Data',sData,'Coords',sCo,'Size',sSize,'Component',sComp);
        
        mapThis = reshape(-1*strain.Data,strain.Size);  %disp('##ISSUE: -1')
                
        switch mapID(1)
            case 1
                mapElements = {mapID(2),1:strain.Size(2),1:strain.Size(3)};
            case 2
                mapElements = {1:strain.Size(1),mapID(2),1:strain.Size(3)};
            case 3
                mapElements = {1:strain.Size(1),1:strain.Size(2),mapID(2)};
        end
        % Coordinates for map
        xx_temp=reshape(strain.Coords(:,x_temp),strain.Size);
        yy_temp=reshape(strain.Coords(:,y_temp),strain.Size);        
        ssX=(squeeze(xx_temp(mapElements{:}))-c0(1));
        ssY=(squeeze(yy_temp(mapElements{:}))-c0(2));        
        
    otherwise
        disp('ERROR: Map value mapID(3) invalid. Valid values are 1-6.')
        ERROR
end

% Map data for required slice
ssMap =(squeeze(mapThis(mapElements{:}))-dispOrigin(mapID(3)));

% Coordinates for vector origins for required slice
vecXData = reshape(data3D_col(:,x_temp+3),dataSize);
vecYData = reshape(data3D_col(:,y_temp+3),dataSize);

sX=(squeeze(xx_temp(vecElements{:}))-c0(1));
sY=(squeeze(yy_temp(vecElements{:}))-c0(2));

% Values for vectors for required slice
sUx=(squeeze(vecXData(vecElements{:}))-d0(1));
sUy=(squeeze(vecYData(vecElements{:}))-d0(2));

%% 2.1 Smooth map

Conv=ssMap;
clear ssMap
ssMap = conv2(Conv,F,'same'); %smoothes map function

%% 2.2 Contour settings
min_col = min(min(ssMap)); max_col = max(max(ssMap)); %Default

% min_col = -6e-3;max_col = 6e-3; disp '##WARNING: Fixed scale'

% min_col = min(min(min(mapThis))); max_col = max(max(max(mapThis))); disp '##WARNING: Fixed scale'

contours = min_col:(max_col-min_col)/CN:max_col;

%% 3. Prepare figure
%figure
clf
hold on

%% 3.1 Plot borders
if PlotBorders
    plotBorder([ssX(:) ssY(:)], mapID)
end

%% 3.2 Plot map
% contourf(ssX',ssY',ssMap',[min_col,contours],'LineStyle','none');

%% 3.3 Plot vectors
if mapID(4)
    h1=quiver(sX,sY,sUx,sUy,V_S);    
    if verLessThan('matlab', '8.4')
        warning('Running on MATLAB version pre 2014.  Upgrade your MATLAB')
        adjust_quiver_arrowhead_size2(h1,A_S);
        set(h1,'linewidth',A_T);
        set(h1,'color',A_C); 
    else
        h1.Color = A_C; %RGB triplet|color string|short color string
        h1.LineWidth = A_T; %in point (1/72")
        h1.MaxHeadSize = A_S;
%         h1.AutoScale = 'off'      %Not generally useful
%         h1.AutoScaleFactor = V_S; %in quiver line
    end
end

%% Define axes and colour bar
axis equal;
%%Choose a pair of suitable axis definitions

% xlim([floor(min(ssX(:,1))) ceil(max(ssX(:,1)))]);
% ylim([floor(min(ssY(1,:))) ceil(max(ssY(1,:)))]);

% xlim([floor(min(ssX(:,1))-0.1*uMult(x_temp)) ceil(max(ssX(:,1))+0.1*uMult(x_temp))]);
% ylim([floor(min(ssY(1,:))-0.1*uMult(y_temp)) ceil(max(ssY(1,:))+0.1*uMult(y_temp))]);

% This pair was to enclose all arrows but does not work properly :(
% xvectors = sX+sUx;yvectors = sY+sUy;
% xlim([floor(min([xvectors(:);sX(:)])-0.1*uMult(x_temp)) ceil(max([xvectors(:);sX(:)])+0.1*uMult(x_temp))]);
% ylim([floor(min([yvectors(:);sY(:)])-0.1*uMult(y_temp)) ceil(max([yvectors(:);sY(:)])+0.1*uMult(y_temp))]);

% xlim([floor(min(ssY(:,1)))-0.5 ceil(max(ssY(:,1)))+0.5]);
% ylim([floor(min(ssX(1,:)))-0.5 ceil(max(ssX(1,:)))+0.5]);

% xlim([floor(min(ssX(1,:)))-0.5 ceil(max(ssX(1,:)))+0.5]);
% ylim([floor(min(ssY(:,1)))-0.5 ceil(max(ssY(:,1)))+0.5]);

% xlim([floor(min(ssX(:,1)+sX(:,1)))-0.5 ceil(max(ssX(:,1)+sX(:,1)))+0.5]);
% ylim([floor(min(ssY(:,1)+sY(1,:)))-0.5 ceil(max(ssY(:,1)+sY(1,:)))+0.5]);

set(gca, 'YDir', 'reverse');
c=colorbar;
colTit = [component{mapID(3)}(1:3) ' ' Title{3} ' ' component{mapID(3)}(4:end)];
ylabel(c,colTit);

%% 3.4 Add labels
title(figTitle);
xlabel(figX);
ylabel(figY);

%% Change font size
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold');

hold off

 end
 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% adjust_quiver_arrowhead_size2
%
function adjust_quiver_arrowhead_size2(quivergroup_handle, scaling_factor)

% Make quiver arrowheads bigger or smaller.
%
% adjust_quiver_arrowhead_size(quivergroup_handle, scaling_factor)
%
% Example:
%   h = quiver(1:100, 1:100, randn(100, 100), randn(100, 100));
%   adjust_quiver_arrowhead_size(h, 1.5);   % Makes all arrowheads 50% bigger.
%
% Inputs:
%   quivergroup_handle      Handle returned by "quiver" command.
%   scaling_factor          Factor by which to shrink/grow arrowheads.
%
% Output: none

% Kevin J. Delaney
% December 21, 2011
% BMT Scientific Marine Services (www.scimar.com)

if ~exist('quivergroup_handle', 'var')
    help(mfilename);
    return
end

if isempty(quivergroup_handle) || any(~ishandle(quivergroup_handle))
    errordlg('Input "quivergroup_handle" is empty or contains invalid handles.', ...
             mfilename);
    return
end

if length(quivergroup_handle) > 1
    errordlg('Expected "quivergroup_handle" to be a single handle.', mfilename);
    return
end

if ~strcmpi(get(quivergroup_handle, 'Type'), 'hggroup')
    errrodlg('Input "quivergroup_handle" is not of type "hggroup".', mfilename);
    return
end

if ~exist('scaling_factor', 'var') || ...
   isempty(scaling_factor) || ...
   ~isnumeric(scaling_factor)
    errordlg('Input "scaling_factor" is missing, empty or non-numeric.', ...
             mfilename);
    return
end

if length(scaling_factor) > 1
    errordlg('Expected "scaling_factor" to be a scalar.', mfilename);
    return
end

if scaling_factor <= 0
    errordlg('"Scaling_factor" should be > 0.', mfilename);
    return
end

line_handles = get(quivergroup_handle, 'Children');

if isempty(line_handles) || (length(line_handles) < 3) || ...
   ~ishandle(line_handles(2)) || ~strcmpi(get(line_handles(2), 'Type'), 'line')
    errordlg('Unable to adjust arrowheads.', mfilename);
    return
end

arrowhead_line = line_handles(2);

XData = get(arrowhead_line, 'XData');
YData = get(arrowhead_line, 'YData');

if isempty(XData) || isempty(YData)
    return
end

%   Break up XData, YData into triplets separated by NaNs.
first_nan_index = find(~isnan(XData), 1, 'first');
last_nan_index  = find(~isnan(XData), 1, 'last');

for index = first_nan_index : 4 : last_nan_index
    these_indices = index + (0:2);
    
    if these_indices(end) > length(XData)
        break
    end
    
    x_triplet = XData(these_indices);
    y_triplet = YData(these_indices);
    
    if any(isnan(x_triplet)) || any(isnan(y_triplet))
        continue
    end
    
    %   First pair.
    delta_x = diff(x_triplet(1:2));
    delta_y = diff(y_triplet(1:2));
    x_triplet(1) = x_triplet(2) - (delta_x * scaling_factor);
    y_triplet(1) = y_triplet(2) - (delta_y * scaling_factor);
        
    %   Second pair.
    delta_x = diff(x_triplet(2:3));
    delta_y = diff(y_triplet(2:3));
    x_triplet(3) = x_triplet(2) + (delta_x * scaling_factor);
    y_triplet(3) = y_triplet(2) + (delta_y * scaling_factor);
    
    XData(these_indices) = x_triplet;
    YData(these_indices) = y_triplet;
end

set(arrowhead_line, 'XData', XData, 'YData', YData);

end

%%%%%%%%%%
%
% sliceVars
% Attempt to tidy up the above function but proved too complex to be useful
%{
function [x_temp,y_temp,vecElements,c0,d0,outOpts,figTitle,figX,figY,DX]...
    = sliceVars(mapID,data3D_col,dataSize,coordinateOrigin,displacementsOrigin,outOpts,Title,dx_mm)

Component = {'Ux';'Uy';'Uz';'Exx';'Eyy';'Ezz'};

if isempty(outOpts)
    %No inputs,so use defaults
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];     %Smoothing function
    CN = 30;                % Number of contours
    %For vector field
    n = 16;                 % Plot every nth vector on a square array
    A_S = 1;                % Arrow head size
    A_C = [1,0,1];          % Vector color
    A_T = 1;                % Arrow thickness
    V_S = 2;                % Arrow size
    outOpts = {F,CN,n,A_S,A_C,A_T,V_S};
end

switch size(dx_mm,2)
    case 0
        disp(' ')
        disp('WARNING: dx_mm not entered - dx_mm = 1 will be used so calculated strain values innaccurate.')
        DX=1;
    case 1
        DX = [dx_mm dx_mm dx_mm];
    case 3
        DX = dx_mm;
    otherwise
        disp('ERROR: Invalid number elements in dx_mm. Valid lengths are 1 or 3')
        ERROR
end

% Define for required slice
switch mapID(1)
    case 1
        x_temp = 2;
        y_temp = 3;
        vecElements = {mapID(2),1:n:dataSize(2),1:n:dataSize(3)};
        figX = 'y /mm';
        figY = 'z /mm';
        Tit1 = 'Y-Z plane (X = ';
    case 2
        x_temp = 1;
        y_temp = 3;
        vecElements = {1:n:dataSize(1),mapID(2),1:n:dataSize(3)};
        figX = 'x /mm';
        figY = 'z /mm';
        Tit1 = 'X-Z plane (Y = ';
    case 3
        x_temp = 1;
        y_temp = 2; 
        vecElements = {1:n:dataSize(1),1:n:dataSize(2),mapID(2)};
        figX = 'x /mm';
        figY = 'y /mm';
        Tit1 = 'X-Y plane (Z = ';
    otherwise
        disp('ERROR: Slice orientation value (mapID(1)) incorrect. Valid values are 1-3.')
        ERROR
end

if Title{1}
    Tit2 = Component{mapID(3)};
    figTitle = [Tit2 ' ' Title{2} ' ' Tit1 num2str(mapID(2)) ') for ' Title{3}];
else
    figTitle = Title{2};
end

switch logical(true)
    case strcmp(coordinateOrigin,'Natural')
        coordOrigin = [0 0 0];
    case strcmp(coordinateOrigin,'Centre')
        coordOrigin = nanmean(data3D_col(:,1:3));
    otherwise
        coordOrigin = coordinateOrigin;
end

switch logical(true)
    case strcmp(displacementsOrigin,'Absolute')
        dispOrigin = zeros(1,6);
    case strcmp(displacementsOrigin,'Relative')
        dispOrigin = [nanmean(data3D_col(:,4:6)) 0 0 0];
    otherwise
        dispOrigin = [displacementsOrigin 0 0 0];
end

c0 = [coordOrigin(x_temp) coordOrigin(y_temp)];
d0 = [dispOrigin(x_temp) dispOrigin(y_temp)];
end
%}
