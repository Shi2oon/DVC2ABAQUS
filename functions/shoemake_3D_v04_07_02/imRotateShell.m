function [dataOut, sizeOut, dxOut] = imRotateShell(data3D_col,dataSize,dx_mm,theta_deg)
%% imRotateShell (data3D_col,dataSize,dx_mm,theta_deg) 
% A shell script to rotate columnwise 3D data. Rotation
% is anti-clockwise about z-axis, centred on image. Also calculates the
% new volume and voxel sizes.
% Inputs:
%   data3D_col      6 column matrix of [coordData displaceData]
%   dataSize        3 element vector of [lx ly lz]
%   dx_mm           Scalar orignal vector size in mm
%   theta_deg       Angle of rotation in degrees
% Outputs:
%   dataOut         As data3D_col
%   sizeOut         As dataSize
%   dxOut           As dx_mm
%
% (c) MSL Jordan, University of Oxford 2015

%%
%{
%TEST VAR
data3D_col = tdata;
% dataSize =    
% dx_mm          
% theta_deg
%}

%% 0. Definitions
if size(theta_deg,2) ==3
    theta_deg = theta_deg(3); 
%     disp '##ISSUE 3-rot/1-rot. Identified 6 Sept 2014';
else
    %do nothing (theta_deg should be 1 value only)
end

%% 1. Prepare columnwise output data matrix (dataOut), and resulting size of 3D volume
temp3D = reshape(data3D_col(:,4),dataSize);
tempRot = imRotCustom(temp3D,theta_deg);

sizeOut = size(tempRot);
dataOut = zeros(numel(tempRot),6);

%% 2. Fill displacement columns of dataOut
dataOut(:,4) = tempRot(:);

for i=5:6
   temp3D = reshape(data3D_col(:,i),dataSize);
   tempRot = imRotCustom(temp3D,theta_deg);
   dataOut(:,i) = tempRot(:);
end

%% 3. Compute new voxel size
switch size(dx_mm,2)
    case 1
        dx_mm3 = repmat(dx_mm,1,3);
    case 3 
%         if dx_mmIn(1) ~= dx_mmIn(2)
%             disp('ERROR: capability of imRotateShell for non-square (prismatic) voxels is not coded. Sorry'); ERROR
%         else
%             dx_mm = dx_mmIn;
%         end
        warning('WARNING: Potential issues with non-regular cubic arrays');
        dx_mm3 = dx_mm;
    otherwise
        error('ERROR: dx_mm is of invalid size. Only 1D and 3D vectors are permitted'); 
end

switch mod(theta_deg, 90) 
    case 0
        dxOut = [dx_mm3(1) dx_mm3(2) dx_mm3(3)];
    otherwise
        dx = dx_mm3(1)/cosd(theta_deg);
        dy = dx_mm3(2)/cosd(theta_deg);
        dz = dx_mm3(3);
        dxOut = [dx dy dz];
        %For none-cubic voxels, calculate dxOut=[dx_mm dy_mm dz_mm]
        %Need to adjust S.4 (see xBinData.m for example)
end

%{
dxD = size(dx_mmIn,2);
switch dxD
    case 1
        dx_mm = dx_mmIn;
        %do nothing
    case 3
        if dx_mmIn(1) ~= dx_mmIn(3)
            disp('ERROR: capability of imRotateShell for non-square (prismatic) voxels is not coded. Sorry'); ERROR
        else
            dx_mm = dx_mmIn(3);
        end
    otherwise
        disp('ERROR: dx_mm is of invalid size'); ERROR
end


switch mod(theta_deg, 90) 
    case 0
        dxOut = dx_mm
    otherwise
        dx = dx_mm/cosd(theta_deg);

        dxOut = [dx dx dx];
        %For none-cubic voxels, calculate dxOut=[dx_mm dy_mm dz_mm]
        %Need to adjust S.4 (see xBinData.m for example)
end
%}
%% 4. Compute new coordinates (first 3 columns of dataOut)

xlin = ((1:sizeOut(1))-0.5)*dxOut(1);
ylin = ((sizeOut(2):-1:1)-0.5)*dxOut(2);
zlin = ((1:sizeOut(3))-0.5)*dxOut(3);
[yy, xx, zz] = meshgrid(ylin,xlin,zlin);

% midX = round((sizeOut(1:2)-1)/2)*dxOut;
% [yy xx zz] = meshgrid(-midX(2):dxOut:midX(2),-midX(1):dxOut:midX(1),(1:sizeOut(3))*dxOut);
dataOut(:,1:3) = [xx(:) yy(:) zz(:)];
%% 5. Create output structure
% 
% outStruct = struct('Data',dataOut,'Size',sizeOut,'DX',dxOut);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% {
%  Function: imRotCustom
%

function varargout = imRotCustom(varargin)
%IMROTATECUSTOM Rotate image (simplified IMROTATE)
%   B = IMROTATECUSTOM(A,ANGLE) rotates image A by ANGLE degrees in a 
%   counterclockwise direction around its center point. To rotate the image
%   clockwise, specify a negative value for ANGLE. IMROTATE makes the output
%   image B large enough to contain the entire rotated image. IMROTATE uses
%   nearest neighbor interpolation, setting the values of pixels in B that 
%   are outside the rotated image to NaN.
%
%   B = IMROTATE(A,ANGLE,METHOD) rotates image A, using the interpolation
%   method specified by METHOD - only 'nearest' is accepted! 
%
%        {'nearest'}  Nearest neighbor interpolation
%
%   B = IMROTATE(A,ANGLE,METHOD,BBOX) rotates image A, where BBOX specifies 
%   the size of the output image B. BBOX is a text string that can have 
%   either of the following values. The default value is enclosed in braces
%   ({}).
%
%        {'loose'}    Make output image B large enough to contain the
%                     entire rotated image. B is generally larger than A.
%
%        'crop'       Make output image B the same size as the input image
%                     A, cropping the rotated image to fit. 
%
%   Class Support
%   -------------
%   The input image can be numeric or logical.  The output image is of the
%   same class as the input image.
%
%   Performance Note
%   ----------------
%   This function may take advantage of hardware optimization for datatypes
%   uint8, uint16, and single to run faster.
%
%   Example
%   -------
%        % This example brings image I into horizontal alignment by
%        % rotating the image by -1 degree.
%        
%        I = fitsread('solarspectra.fts');
%        I = mat2gray(I);
%        J = imrotate(I,-1,'bilinear','crop');
%        figure, imshow(I), figure, imshow(J)
%
%   See also IMCROP, IMRESIZE, IMTRANSFORM, TFORMARRAY.

%   Copyright 1992-2009 The MathWorks, Inc.
%   $Revision: 5.25.4.9 $  $Date: 2009/05/14 16:58:13 $
%   $Edited by M Jordan 2014$ 

% Grandfathered:
%   Without output arguments, IMROTATE(...) displays the rotated
%   image in the current axis.  

[A,ang,method,bbox] = parse_inputs(varargin{:});

so = size(A);
twod_size = so(1:2);

if mod(ang,360) ==0;
  % Catch zero rotation
      B = A;

else % Perform general rotation
    
    phi = -ang*pi/180; % Convert to radians
    
    rotate = maketform('affine',[ cos(phi)  sin(phi)  0; ...
        -sin(phi)  cos(phi)  0; ...
        0       0       1 ]);
    
    [loA,hiA,loB,hiB,outputSize] = getOutputBound(rotate,twod_size,bbox);
    
    % rotate using tformarray
                
    boxA = maketform('box',twod_size,loA,hiA);
    boxB = maketform('box',outputSize,loB,hiB);
    T = maketform('composite',[fliptform(boxB),rotate,boxA]);

    R = makeresampler('linear','fill');
    B = tformarray(A, T, R, [1 2], [1 2], outputSize, [], nan);       
end
   
% Output
switch nargout,
    case 0,
      wid = 'Images:imrotate:obsoleteSyntax';    
      warning(wid, '%s', ['Obsolete syntax. In future versions IMROTATE ',... 
      'will return the result in ans instead of displaying it in figure.']);
      imshow(B);
    case 1,
      varargout{1} = B;
    case 3,
      wid = 'Images:imrotate:obsoleteSyntax';    
      warning(wid, '%s', ['[R,G,B] = IMROTATE(RGB) is an obsolete output syntax. ',...
      'Use one output argument to receive the 3-D output RGB image.']);
          for k=1:3,
            varargout{k} = B(:,:,k);
          end;
    otherwise,
      eid = 'Images:imrotate:tooManyOutputs';    
      error(eid, '%s', 'Invalid number of output arguments.');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function: getOutputBound
%
function [loA,hiA,loB,hiB,outputSize] = getOutputBound(rotate,twod_size,bbox)

% Coordinates from center of A
hiA = (twod_size-1)/2;
loA = -hiA;
if strcmpi(bbox, 'loose')  % Determine limits for rotated image
    hiB = ceil(max(abs(tformfwd([loA(1) hiA(2); hiA(1) hiA(2)],rotate)))/2)*2;
    loB = -hiB;
    outputSize = hiB - loB + 1;
else % Cropped image
    hiB = hiA;
    loB = loA;
    outputSize = twod_size;
end
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function: parse_inputs
%

function [A,ang,method,bbox] = parse_inputs(varargin)
% Outputs:  A       the input image
%           ang     the angle by which to rotate the input image
%           method  interpolation method (nearest,bilinear,bicubic)
%           bbox    bounding box option 'loose' or 'crop'

% Defaults:
method = 'n';
bbox = 'l';

narginchk(2,4);
switch nargin
case 2,             % imrotate(A,ang)        
  A = varargin{1};
  ang=varargin{2};
case 3,             % imrotate(A,ang,method) or
  A = varargin{1};  % imrotate(A,ang,box)
  ang=varargin{2};
  method=varargin{3};
case 4,             % imrotate(A,ang,method,box) 
  A = varargin{1};
  ang=varargin{2};
  method=varargin{3};
  bbox=varargin{4};
otherwise,
  eid = 'Images:imrotate:invalidInputs';    
  error(eid, '%s', 'Invalid input arguments.');
end

% Check validity of the input parameters 
if ischar(method) && ischar(bbox),
  strings = {'nearest','bilinear','bicubic','crop','loose'};
  idx = strmatch(lower(method),strings);
  if isempty(idx),
    eid = 'Images:imrotate:unrecognizedInterpolationMethod';
    error(eid, 'Unknown interpolation method: %s', method);
  elseif length(idx)>1,
    eid = 'Images:imrotate:ambiguousInterpolationMethod';
    error(eid, 'Ambiguous interpolation method: %s', method);
  else
    if idx==4,bbox=strings{4};method=strings{1};
    elseif idx==5,bbox = strings{5};method=strings{1};
    else method = strings{idx};
    end
  end  
  idx = strmatch(lower(bbox),strings(4:5));
  if isempty(idx),
    eid = 'Images:imrotate:unrecognizedBBox';
    error(eid, 'Unknown BBOX parameter: %s', bbox);
  elseif length(idx)>1,
    eid = 'Images:imrotate:ambiguousBBox';
    error(eid, 'Ambiguous BBOX string: %s', bbox);
  else
    bbox = strings{3+idx};
  end 
else
  eid = 'Images:imrotate:expectedString';
  error(eid, '%s', 'Interpolation method and BBOX have to be a string.');  
end
end
%}
