function [data_array, wrnmsg] = reformData(data_col, dataSize)
% Reforms columnwise 3D data (xyz) into volumewise 4D matrix. Assumes
% coordinates are created as:
%
% [yy,xx,zz] = MESHGRID(ymax:-1:ymin,xmin:xmax,zmin:zmax)
%
% See also reformColumns
%
% (c) 2015, MSL Jordan, University of Oxford and EDF

% disp '##WARNING: reformData modified 1 Sept 2014 @ 12.28. Check output'
%% 0. Definition
cols = size(data_col,2);
data_array = zeros([dataSize,cols]);

%% 1. Reshape data
switch cols
    case {2,4}
        for q = 1:cols
            data_array(:,:,q) = reshape(data_col(:,q),dataSize);    %3D
        end
        wrnmsg = 'data_array is 3D as 2D data was detected';
    case {3,6}
        for q = 1:cols
            data_array(:,:,:,q) = reshape(data_col(:,q),dataSize);  %4D
        end
        wrnmsg = '';
    otherwise
        for q = 1:cols
            data_array(:,:,:,q) = reshape(data_col(:,q),dataSize);  %4D
        end
        wrnmsg = 'N-Dimensioanl data detected; data_array is 4D';
end
