function [ corr_fea_data ] = AdjustmentBr0_v1beta(OriPts,DefPts,Def_data)
%Selim Barhli 2014
%Adjust a 2D or 3D  FE dataset in the DVC axis system.
%use : [ corr_fea_data ] = AdjustmentBr0_v1beta(OriPts,DefPts,Def_data)
%
%OriPts is an array of 2/3 points (x,y,z) in the original image (DIC/DVC)
%DefPts is the coordinates of those 2/3 points in the image to adjust (FE
%model).
%Def_data contains in a 4/6 columns data style the displacement field of the
%FE model.

% fprintf('[*] AdjustmentBr0 plugin : ON.\n');
fprintf('[*] AdjustmentBr0 plugin : 935N.\n');


%get the dimentionality of the dataset
dime = length(OriPts);
if (length(DefPts) ~= dime || size(Def_data,2)~= 2*dime)
    fprintf('[ERROR] in AdjustmentBr0 - Please check dimensionality consistency.\n')
end

%calculate vector from points
if(dime == 2)
    vecOri = diff(OriPts);
    vecDef = diff(DefPts);
elseif (dime==3)
    vecOri = diff(OriPts(1:2,:));
    vecDef = diff(DefPts(1:2,:));
    vecOri_2 = diff(OriPts([1 3],:));
    vecDef_2 = diff(DefPts([1 3],:));
    
end


%%Correcting for scaling issues
%determining norm of the ref vectors in the two images
fprintf('\t[*] Correcting FE model for scaling : ');

normOri = norm(vecOri); normDef = norm(vecDef);
normOri_2 = norm(vecOri_2);normDef_2 = norm(vecDef_2);

%calculating the scaling factor
scale_f = normOri/normDef;

%transform the FEdata to have the same scale than DIC data and updating the
%ref vector in the FE data;
Def_data = Def_data.*scale_f;
DefPts = DefPts.*scale_f;
vecDef = vecDef.*scale_f; vecDef_2 = vecDef_2.*scale_f;
normDef = norm(vecDef);
fprintf('OK\n');


%%Correcting for rotation issues

%in the 2D case, getting the angle between the two vector is sufficient. In
%3D, first correct the first vector for rool, pitch and yaw. Then compute
%the solid angle between the second vectosr in the plane where the first
%vectors are the normal. Finally corect using the rotation matrix "from
%axis and angle".

fprintf('\t[*] Correcting FE model for rotation : ');

if (dime == 2)
    theta = -(atan2(vecDef(2),vecDef(1)) - atan2(vecOri(2),vecOri(1)));
    matrot = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    
    %correct the rotation and update the ref vector
    Def_data(:,1:2)=(matrot*Def_data(:,1:2)')';
    Def_data(:,3:4)=(matrot*Def_data(:,3:4)')';
    DefPts = (matrot*DefPts')';
    vecDef = diff(DefPts); normDef = norm(vecDef);
elseif (dime == 3)
    %correct separatly vect1 for yaw, pitch and roll
    
    %roll correction
    rDef = vecDef(2:3); rOri = vecOri(2:3);
    if (sum(rDef==0)==2 | sum(rOri==0)==2)
        thetaroll = 0;
    else
        thetaroll = -(atan2(rDef(2),rDef(1)) - atan2(rOri(2),rOri(1)));
        %thetaroll(thetaroll==pi())=0;
    end
    matroll = [1 0 0;0 cos(thetaroll) -sin(thetaroll); 0 sin(thetaroll) cos(thetaroll)];
    vecDef = (matroll*vecDef')';
    
    
    %pitch correction
    pDef = vecDef([1 3]); pOri = vecOri([1 3]);
    if (sum(pDef==0)==2 | sum(pOri==0)==2)
        thetapitch = 0;
    else
        %VERIFIE
        thetapitch =(atan2(pDef(2),pDef(1)) - atan2(pOri(2),pOri(1)));
        %thetapitch(thetapitch==pi())=0;
    end
    matpitch = [cos(thetapitch) 0 sin(thetapitch);0 1 0;-sin(thetapitch) 0 cos(thetapitch)];
    vecDef = (matpitch*vecDef')';
    
    %yaw correction
    yDef = vecDef(1:2); yOri = vecOri(1:2);
    if (sum(yDef==0)==2 | sum(yOri==0)==2)
        thetayaw = 0;
    else
        %VERIFIE
        thetayaw = -(atan2(yDef(2),yDef(1)) - atan2(yOri(2),yOri(1)));
        %thetayaw(thetayaw==pi())=0;
    end
    matyaw = [cos(thetayaw) -sin(thetayaw) 0;sin(thetayaw) cos(thetayaw) 0;0 0 1];
    vecDef = (matyaw*vecDef')';
    
    %MODULO TO CHECK FOR 180
%     thetaroll = mod(thetaroll,2*pi());
%     thetapitch = mod(thetapitch,2*pi());
%     thetayaw = mod(thetayaw,2*pi());
    
    %matroll = [1 0 0;0 cos(thetaroll) -sin(thetaroll); 0 sin(thetaroll) cos(thetaroll)];
%     matpitch = [cos(thetapitch) 0 sin(thetapitch);0 1 0;-sin(thetapitch) 0 cos(thetapitch)];
%     matyaw = [cos(thetayaw) -sin(thetayaw) 0;sin(thetayaw) cos(thetayaw) 0;0 0 1];
    
    matrot = matyaw*matpitch*matroll;
    
    %correct the rotation for vect1 and update the ref vectors
    Def_data(:,1:3)=(matrot*Def_data(:,1:3)')';
    Def_data(:,4:6)=(matrot*Def_data(:,4:6)')';
    DefPts = (matrot*DefPts')';
    
    vecDef = diff(DefPts(1:2,:));
    vecDef_2 = diff(DefPts([1 3],:));
    
    %Now find the oriented solid angle between vect2 in the plane where vect1 is
    %normal
    unitplan = vecDef./norm(vecDef);
    projDef2 = cross(unitplan,cross(vecDef_2,unitplan));
    projOri2 = cross(unitplan,cross(vecOri_2,unitplan));
    
    solangl = acos(dot(projDef2,projOri2)/(norm(projDef2)*norm(projOri2)));
    orientsgn = sign(det([projDef2;projOri2;unitplan]));
    solangl = real(solangl)*orientsgn;
%     solangl(solangl==pi())=0;
    
    crossvecDef = [0 -unitplan(3) unitplan(2);unitplan(3) 0 -unitplan(1);-unitplan(2) unitplan(1) 0];
    
    matrot2 = cos(solangl)*eye(3)+sin(solangl)*crossvecDef+(1-cos(solangl))*kron(unitplan,unitplan');
    
    Def_data(:,1:3)=(matrot2*Def_data(:,1:3)')';
    Def_data(:,4:6)=(matrot2*Def_data(:,4:6)')';
    DefPts = (matrot2*DefPts')';
    
    vecDef = diff(DefPts(1:2,:));
    vecDef_2 = diff(DefPts([1 3],:));
end
fprintf('OK\n');

%%Correcting for RBM
fprintf('\t[*] Correcting FE model for RBM : ');
rbmvec = median(DefPts-OriPts);
% rbmvec = 
% rbmvec = roundn(DefPts(1,:)-OriPts(1,:),-4)
DefPts = DefPts-repmat(rbmvec,size(DefPts,1),1);
Def_data(:,1:dime) = Def_data(:,1:dime)- repmat(rbmvec,size(Def_data,1),1);
fprintf('OK\n');

corr_fea_data = Def_data;
