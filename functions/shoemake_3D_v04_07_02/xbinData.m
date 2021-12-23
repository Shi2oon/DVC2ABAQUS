function [dataOut, sizeOut, dxOut] = xbinData(request,data_col,dataSize,dx_mm,nBins)    
% (c) MSL Jordan, University of Oxford 2015

%% Process
%{
0. Definitions
1. Reshape data and bin
1. Mirror data in {x-z} plane  - answers on a postcard please!
2. Calculate dx and create new coordinate set
3. Form dataOut 
4. Calc slice-wise average displacements
%}

%% Definitions
switch size(data_col,2)
    case 4
        dims = 2;
        %Generate 3D dataset: Inefficient method of working with 2D data
        datLen = size(data_col,1);
        col3 = [ones(datLen,1); 2*ones(datLen,1)];
        col6 = zeros(2*datLen,1);
        compData = repmat(data_col,2,1);

        data3D_col = [compData(:,1:2) col3 compData(:,3:4) col6];
        dataSize = [dataSize(1:2) 2];
    case 6
        dims = 3;
        data3D_col = data_col;
    otherwise
        disp '##ERROR: data_col of invalid size. Valid number of columns is 4 (2D) or 6(3D)##';ERROR
end

binSize = floor(dataSize(1)/nBins);
sizeOut = [nBins dataSize(2) dataSize(3)];

crUx_xbin = zeros(sizeOut);
crUy_xbin = zeros(sizeOut);
crUz_xbin = zeros(sizeOut);

if request ~= 1
    disp('##ERROR:I do not like that direction. Try again :p ##'); ERROR
end

switch size(dx_mm,2)
    case 0
        disp(' ')
        disp('WARNING: dx_mm not entered - dx_mm = 1 will be used so calculated bin size and strain values innaccurate.')
        DX=[1 1 1];
    case 1
        DX = [dx_mm dx_mm dx_mm];
    case 3
        DX = dx_mm;
    otherwise
        disp('ERROR: Invalid number elements in dx_mm. Valid lengths are 1 or 3'); ERROR
end

%% 1. Reshape and bin data

Ux3D = reshape(data3D_col(:,4),dataSize);
Uy3D = reshape(data3D_col(:,5),dataSize);
Uz3D = reshape(data3D_col(:,6),dataSize);

for i = 1:nBins
    fSlice = (i-1)*binSize+1;
    lSlice = i*binSize;
    crUx_xbin(i,:,:) = nanmean(Ux3D(fSlice:lSlice,:,:),1);
    crUy_xbin(i,:,:) = nanmean(Uy3D(fSlice:lSlice,:,:),1);
    crUz_xbin(i,:,:) = nanmean(Uz3D(fSlice:lSlice,:,:),1);
end

%% 2. Calc dxOut and new coordinates

dxOut = [DX(1)*binSize DX(2) DX(3)];

xlin = ((1:sizeOut(1))-0.5)*dxOut(1);
ylin = ((sizeOut(2):-1:1)-0.5)*dxOut(2);
zlin = ((1:sizeOut(3))-0.5)*dxOut(3);
[ybin, xbin, zbin] = meshgrid(ylin,xlin,zlin);

%% 3. Form dataOut

dataOut = [xbin(:) ybin(:) zbin(:) crUx_xbin(:) crUy_xbin(:) crUz_xbin(:)];

if dims ==2
    if numel(unique(dataOut(:,3)))~=2
        disp '##ERROR: The permutations are screwed-up. Contact MSLJ##'; ERROR
    end
    dataOut = dataOut(1:datLen,:);
    dataOut(:,3:3:6) = [];
end

%% 4. Calc slice-wise average displacements
% %     Ux2DAve = zeros(numberOfBins,1);
%     Ux2DAve = nanmean(nanmean(cleaned_Uy_rot(fSlice:lSlice,:,:)));
%     plot(squeeze(Ux2DAve)) 
    
end