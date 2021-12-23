function [dataOut,coordsOut,sizeOut,compOut] = calcOrthoStrain(strainID,data3D_col,dataSize,method,dx_mm)
% (c) MSL Jordan, University of Oxford 2015

%% Process
%{
0. Definitions
1. Create Dummy 3D matirices of required displacement
2. Average dummy matrices perpendiculat to orthogonal direction
3. Subtract displacments along required direction to extract relative displacements
4. Divide by voxel spacing to calc strain
%}
%% 0. Definitions

if size(data3D_col,2) == 4
    dims = 2;
    %Generate 3D dataset: Inefficient method of working with 2D data
    datLen = size(data3D_col,1);
    col3 = [ones(datLen,1); 2*ones(datLen,1)];
    col6 = zeros(2*datLen,1);
    compData = repmat(data3D_col,2,1);
    
    data3D_col = [compData(:,1:2) col3 compData(:,3:4) col6];
    dataSize = [dataSize(1:2) 2];
else 
    dims =3;
end

switch size(dx_mm,2)
    case 0
        disp(' ')
        disp('WARNING: dx_mm not entered - dx_mm = 1 will be used so calculated strain values innaccurate.')
        DX=[1 1 1];
    case 1
        DX = [dx_mm dx_mm dx_mm];
    case 3
        DX = dx_mm;
    otherwise
        disp('ERROR: Invalid number elements in dx_mm. Valid lengths are 1 or 3'); ERROR
end

%% 1. Permute data

switch strainID
    case 4
        tempU = reshape(data3D_col(:,strainID),dataSize);
        tempU = permute(tempU,[2 3 1]);
        tempDX = permute(DX,[2 3 1]);
        compOut = 'Exx';
        
    case 5
        tempU = reshape(data3D_col(:,strainID),dataSize);
        tempU = permute(tempU,[3 1 2]);
        tempDX = permute(DX,[3 1 2]);
        compOut = 'Eyy';
    case 6
        tempU = reshape(data3D_col(:,strainID),dataSize);
        tempDX = DX;
        compOut = 'Ezz';
        if dims == 2
            disp '##WARNING: Ezz is undefined with 2D data##'
        end
    otherwise
        disp('ERROR Invalid strainID. Valid values are 4-6 for Exx-Ezz respectively'); ERROR
end

tempSize = size(tempU);
lx = tempSize(1);
ly = tempSize(2);
lz = tempSize(3);

%% 2. Dummy matrix strain calculation
switch method
    case 'Simple'
        dummy1=zeros(lx,ly,lz+1);
        dummy2=zeros(lx,ly,lz+1);
        dummy1(:,:,1:lz)=tempU;
        dummy2(:,:,2:lz+1)=tempU;
        tempEz=(dummy2-dummy1)./tempDX(strainID-3);
        tempOut=tempEz(:,:,1:lz);
        
        sizeOut = dataSize;
        coordsOut = data3D_col(:,1:3);        
    case '2x2x2'
        for i= 1:2
            dummy1=zeros(lx,ly,lz+1);
            dummy2=zeros(lx,ly,lz+1);
            dummy3=zeros(lx,ly,lz+1);
            dummy4=zeros(lx,ly,lz+1);
            % dummy3
            % dummy4
            % dummy5
            % dummy6
            % dummy7
            % dummy8
            avDumSlice = (dummy1+dummy2+dummy3+dummy4)/4;
        end
        disp('##ERROR: Options other than "Simple" not coded. Sorry')
        ERROR
        tempOut=tempEz(:,:,1:lz);
        sizeOut = dataSize;
        coordsOut = data3D_col(:,1:3);
    otherwise
        disp('##ERROR: Options other than simple Uz not coded. Sorry')
        ERROR
end

%% 3. Un-permute data
switch strainID
    case 4
        tempOut = permute(tempOut,[3 1 2]);
    case 5
        tempOut = permute(tempOut,[2 3 1]);
    case 6
        %Nothing happens
end

dataOut = tempOut(:);

%% 4. Collapse 2D data
if dims ==2
    if numel(unique(dataOut(:,3)))~=2
        disp '##ERROR: The permutations are screwed-up. Contact MSLJ##'; ERROR
    end
    
    dataOut = dataOut(1:datLen,:);
    dataOut(:,3:3:6) = [];
end

%% Create struct to keep data together
% 
% strain_struct = struct('Data',sData(:),'Coords',dummyCoords,'Size',dummySize,'Component',comp);
end