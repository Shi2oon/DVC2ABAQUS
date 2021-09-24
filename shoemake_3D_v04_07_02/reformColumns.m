function [data3D_col,dataSize] = reformColumns(data_4D)
% Reforms volumewise 4D matrix in columnwise 3D data (xyz). Assumes coordinates are created
% as:
%
% [yy,xx,zz] = MESHGRID(ymax:-1:ymin,xmin:xmax,zmin:zmax)
%
% (c) 2014, MSL Jordan, University of Oxford and EDF

%% 0. Definitions
numPts = numel(data_4D(:,:,:,1));
data3D_col = zeros(numPts,size(data_4D,4));
%% 1. Data size
dataSize = size(squeeze(data_4D(:,:,:,1)));
%% 2. Coordinates
for i = 1:size(data_4D,4)
    data3D_col(:,i) = reshape(data_4D(:,:,:,i),numPts,1);
end
end