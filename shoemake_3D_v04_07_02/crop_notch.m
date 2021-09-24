function [tdata_mask,c_left,theta_deg, symPlanes, alignPts, c0Req] = crop_notch(tdata,notchDiameter,dataSize,dx_pix)
% crop_notch masks the volumetric dataset, in a predefined geometric
% pattern.
%
% Written for EE8519 experiments, the function also outputs veriables
% related to that analysis
%
% (c) MSL Jordan, University of Oxford 2015

%% Process
% 1. generate ones matrix for (originally: top 30) all slices and delete elements in notch = mask [Size Notch_diameter]
% 2. Reshape mask = mask_col [Size]
% 3. Muliply = tdata_mask [tdata]
%}
%% 1.
% smallSize = dataSize;        %[dataSize(1:2) dataSize(3)];
switch notchDiameter;
    case 2
    disp '##ISSUE: Check 2mm notch location';
        c_left = [1 794 -740]/dx_pix; %makeOval(194, -955, 1232, 1232);
        c_left = ([1 754 -740]-1)/dx_pix+1; %makeOval(194, -955, 1232, 1232);
        dec_deg = 0;     %actual value -1.37, but using 0 due to non-cylindrical notch;
        theta_deg = 6.5; %6.25;
        r = 1000/dx_pix;    
    case 4
        c_left = [1 88 -146];   %for dx = 12pixel
        c_left = (c_left-1)/dx_pix*12 +1;
        dec_deg = -0.514;
        theta_deg = 8.1; %9.15
        r = 2000/dx_pix;      

    otherwise
        disp 'WARNING: Invalid Notch_diameter. Valid values are 2 or 4.'
        disp 'Notch has not been cropped.' 
        disp 'Press any key to continue...'
        pause
end

mask3D = notchMask(dataSize,c_left,dec_deg,theta_deg,r);

%% 3. 
L = prod(dataSize);
mask_col = reshape(mask3D,L,1);
mask = [mask_col mask_col mask_col];

%% 4.
tdata_mask = tdata;
tdata_mask(1:L,4:6) = tdata_mask(1:L,4:6).*mask;

%% Check output
%{
x_temp = reshape(tdata_mask(:,4),Size);
figure
M = squeeze(x_temp(:,:,10)); contourf(M); axis image, shg
%}

%% 5. Prepare extra data for other functions
switch notchDiameter
    case 2
        %regPoints
        % Alignment points (3 row vectors of coordinates)
        DVC_points = [...
            0.3024,3.141 ,0.54 ;... %currently for EE8519 2mm
            5.0184,2.6316,0.4824 ;...
            5.0184,2.6316,4.4824];

        FEM_points = [...                %currently for EE8519 2mm
            -2.5,0,12;...          
             2.5,0,12;...
             2.5,0,16];

%         c0Req = [2.6604 3.4623 0.5112];
%         annotOpt                      %add cirles to mark edge of sample
    case 4
        %regPoints
        % Alignment points (3 row vectors of coordinates)
        DVC_points = [...
            0.5922 ,2.7306 ,0.4770 ;... %currently for EE8519 4mm
            5.2812,1.98,0.3924 ;...
            5.2956,1.9782,2.6478];

        FEM_points = [...                %currently for EE8519 4mm
            -2.5,0,12;...          
             2.5,0,12;...
             2.5,0,14.3755];

%          c0Req = [2.9367 2.3553 0.4347];  %using mean of registration
%          points
%         c0Req = [3.121 2.5666 0.288];     %measured
%         annotOpt = 
end
alignPts = [DVC_points; FEM_points];

plane1 = [0,0,0;... %currently y-z plane
          0,1,0;...
          0,0,1];
plane2 = []; 
plane3 = [];
symPlanes = [plane1; plane2; plane3];

c0Req = mean(DVC_points(1:2,:)); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mask = notchMask(dataSize,c_left,dec_deg,theta_deg,r)
%% notch_mask (alpha)
% 
% Defines a 3-D cylindical masking function specifically for masking blunt notches
% cut in the upper x-y planes of a sample, along an arbitrary central axis.
%
% Size = [lx ly lz] : the required mask size
% c_left =[yc zc] : the centre of the circle projected on the slice X=1
% dec = declination angle (degrees): the angle the notch makes with the horizontal
% in a rotated slice, starting from X'=1
% az  =azimuthal angle (degrees): the notch angle against the x = 0 in a Z=1
% slice, starting from X'=1
% r = radius : notch radius as projected on a Z=1 slice. If r=0, calculated
% from Notch_diameter

%{ 
%TEST VAR
Notch_diameter = 4;
Size = [267 213 100];
c_left = [1 61 -146]; %x-z coordinate of centre on Y = 1 plane unrotated plane
dec = -0.514;        %declination notch angle (degrees) in y-z rotated plane
az = 8.1;          %azimuthal notch angle (degrees) in x-y unrotated plane

%}
%% process
%{
% for x-z slice
% 1. Calculate centre = c_temp [c_right, dec, azimuth]
% 2. calc corrected x_row and z_col = x_row, z_col [Size, c_temp] 
% 3. combine x_row and z_col quadratically = radii_matrix [x_row z_col]
% end
% 4. logic (< r ==0) = slice_mask [radii_matrix]
%}

%% Definitions
az = deg2rad(theta_deg);
dec = deg2rad(dec_deg);
% if r==0           %this is almost entirely rubbish
%     r = notchDiameter/2 / (12*1.8e-3) *cos(az);   %projection of notchRadius * 46.3 vox per mm on x-z plane
% end
radii_matrix = zeros(dataSize);
% c_temp = zeros(Size(2),2);

 for j = 1:dataSize(1)  %y-z slices %voxel
     %% 1.
     c_temp = [(c_left(2) + j*tan(az)) (c_left(3) + j*tan(dec))];
     
     %% 2.
     y_row = (1:dataSize(2))- c_temp(1);
     z_col = (1:dataSize(3))- c_temp(2);
     %% 3.
     yy_temp = repmat(y_row.^2,dataSize(3), 1);
     zz_temp = repmat(z_col.^2,dataSize(2), 1);
     radii_matrix(j,:,:) = yy_temp'+zz_temp;
 end
 
%% 4.
mask = round(radii_matrix);
mask(mask <= r^2)=0;
mask(mask > r^2)=1;   %This is now a masking function
end