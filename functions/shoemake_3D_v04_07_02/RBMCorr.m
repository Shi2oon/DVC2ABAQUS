function [dataOut,eulerAngles,rotCentre,RBD] = RBMCorr(centre,data3D_col,dataSize,dx_mm)
% RBMCORR removes the 3D rigid body motions (translation and rotation from 
%     a volumetric dataset (e.g. from DVC). The algorithm for rotation 
%     removal is based on Ken Shoemake's Euler angle correction.
%
% Inputs:
%     dat3D_col       6 column data as [coordinates displacements], as from
%                     DaVis - originally 'cdata'
%     dataSize        3 element vector giving size of dataset [X Y Z]
%     
% Outputs:
%     dataOut         6 column data corrected for rigid body movement
%     eulerAngles     calculated Euler angles 
%     rotCentre       calculated rotation centre
%     
% Code originally developed by M. Mostafavi, with additions by M. Jordan.
% Based on Ken Shoemake's Euler angle extraction.
%
% Last edit: November 2014 (M. Jordan)
% (C) 2014, University of Oxford

disp(' ')
%% 0. Input data check and definitions
%not implemented

%% 1. Clean the dataset

cdata=data3D_col(~any(isnan(data3D_col),2),:);

if(sum(sum(sum(isnan(cdata))))) ~=0
    disp('WARNING: data contains vectors with partial NaN entries');
end

%% 2. Subtract RBD from displacements
C0=zeros(1,6);
RBD = nanmean(cdata(:,4:6));      %row of 3 RBD coordinates
C0(1,4:6)= RBD;
fprintf('RBD /mm: (% 1.4f,% 1.4f,% 1.4f).\n',RBD)
cdata = bsxfun(@minus,cdata, C0);
%bsxfun performs same operation as following two lines
% U0=repmat(C0,size(cdata,1),1);    
% cdata=cdata-U0;

%% 3. Check for rotation correction request
if ~strcmpi(centre,'norot')
    %% 4. Finding the centre of rotation and setting as coordinate origin
    % Use one of (2) methods to calculate position of rigid body motion and set
    % as coordinate origin.
    %
    % X0 coords = index coords of rotn centre

    switch centre
        case 'Minimum'
            %For using centre as minimum displacement vector
            [~,X0L]=min(sum(abs(cdata(:,4:6)')));     %minimum vector of cdata
    %         yy = sum(abs(cdata(:,4:6)'));testsize = size(yy,2);xx=1:testsize;
    %         p = polyfit(xx,yy,2);fit = polyval(2,xx);figure; plot(xx,yy,'.',xx,fit,'g-',X0L,v,'rd');
            rotCentre=cdata(X0L,1:3);                 %X0 coords = index coords of rotn centre
%             disp('##ISSUE: Improve minimum finding algorithm##')
        case 'True'
            %For using centre after RBD
            rotCentre = (nanmean(cdata(~isnan(sum(cdata(:,4:6),2)),1:3)));
        case 'Boundary centre'
            rotCentre=ceil(dataSize./2).*dx_mm;    %X0 coords = index coords of rotn centre
        otherwise
            fprintf('##ERROR: Invalid centre request string. \n###Valid requests are: \n\t"Minimum"\n\t"Boundary centre"\n\t"True"'); ERROR
    end

    %Set rotn centre as coord origin
    X0 = [rotCentre 0 0 0];
    cdata = bsxfun(@minus,cdata,X0);

    %xi0 is the location of points before loading and xip is the location after
    %displacement

    xi0=cdata(:,1:3);                       % Reference coords
    xip=cdata(:,1:3)+cdata(:,4:6);      % Comparison dataset coords

    %% 4.1 Calculating rotation matrix, and checking the rank 
    % Rank should equal 3 to avoid singularity
    rank_xi0 = rank(xi0');
    rank_xip = rank(xip');

    if rank_xi0 ~= 3 || rank_xip ~= 3
        error('xi0 or xip rank incorrect. rank_xi0 = %i rank_xip = %i',rank_xi0,rank_xip);
    end

    R = xi0\xip;  % backslash is MATLAB matrix division - solves for rotation efficiently

    %% 4.2 Extracting Euler angles from the rotation matrix
    %{
    % Theoretical rotation matrix based on the calculated angles: 
    t = [psi theta phi]
    c = cos(t); s = sin(t);
    R_theo=   [c(2)c(3), s(1)s(2)c(3)-c(1)s(3), c(1)s(2)c(3)+s(1)s(3);...
               c(2)s(3), s(1)s(2)s(3)+c(1)c(3), c(1)s(2)s(3)-s(1)c(3);...
               -s(2),    s(1)c(2),              c(1)c(2)            ];
    N.B. There are always at least 2 solutions to t.
    For more info see http://tinyurl.com/kg3ehzq
    %}
    if abs(R(3,1))~=1
        theta = -asin(R(3,1));
    %     theta2 = pi - theta;   %This is the redundant 2nd solution;
    %                            %subtitution of theta generates the 3
    %                            %angles. The use of both angles may aid
    %                            %refinement (not implemented)
        psi = atan2(R(3,2)/cos(theta),R(3,3)/cos(theta));
        phi = atan2(R(2,1)/cos(theta),R(1,1)/cos(theta));
    else
        phi = 0;
        if R(3,1) == -1
            theta = pi/2;
            psi = phi+atan2(R(1,2),R(1,3));
        else
            theta = -pi/2;
            psi = -phi+atan2(-R(1,2),-R(1,3));
        end
    end

    euler = [psi theta phi];

    % disp(['euler = ' num2str(euler,'% 1.3f ')])
    % disp(['Euler angle estimate /degree: (gamma) = (' num2str(radtodeg(euler),'% 1.3f ') ').'])

    %% 4.3.1 Initial estimate of Rotation matrix 
    %OC2
    %% 4.3.2 Refine estimate of Euler angles and Rot matrix
    disp '>>Solver details (refining euler angles):'
    disp '***********'

    %Check Matlab version as optimoptions.m released with 2013a
    MATVer = version('-release');
    if str2double(MATVer(1:4))>= 2013
        options2 = optimoptions('lsqnonlin','MaxFunEval',1E9,'MaxIter',1E12,'TolFun',1e-10,'Display','final'); %'Algorithm','sqp'
    else
        warning('Update MATLAB version to 2014 or above. Using default settings for solver');
        optPath = 'RBMCorr_solverOptions.mat';
        if exist(optPath,'file')
            defaultOpts = load(optPath);
            options2 = defaultOpts.options2;
        else
            error('##ERROR: Default options setting file <%s> not found',optPath); 
        end    
    end

    s = @(x) sin(x);
    c=@(x) cos(x);
    thetaSolver3D = @(t) (R - ...
        [c(t(2))*c(t(3)), s(t(1))*s(t(2))*c(t(3))-c(t(1))*s(t(3)), c(t(1))*s(t(2))*c(t(3))+s(t(1))*s(t(3));...
         c(t(2))*s(t(3)), s(t(1))*s(t(2))*s(t(3))+c(t(1))*c(t(3)), c(t(1))*s(t(2))*s(t(3))-s(t(1))*c(t(3));...
         -s(t(2)),        s(t(1))*c(t(2)),                         c(t(1))*c(t(2))              ]);

    [euler,~] = lsqnonlin(thetaSolver3D,euler,[-2*pi;-2*pi;-2*pi], [2*pi;2*pi;2*pi],options2);
    disp '***********'
    disp '>>End Solver details'
    disp ' '

    c=cos(euler); s=sin(euler);
    R_sol=[c(2)*c(3), s(1)*s(2)*c(3)-c(1)*s(3), c(1)*s(2)*c(3)+s(1)*s(3);...
            c(2)*s(3), s(1)*s(2)*s(3)+c(1)*c(3), c(1)*s(2)*s(3)-s(1)*c(3);...
            -s(2),     s(1)*c(2),                c(1)*c(2)              ];

    % disp(['euler = ' num2str(euler,'% 1.3f ')])
    disp(['Euler angle /degree: (alpha beta gamma) = (' num2str(radtodeg(euler),'% 1.3f ') ').'])

    %{
    disp(['det(R_inital) = ' num2str(det(R))])
    disp('R_initial = ') 
    disp(num2str(R))
    disp(' ')
    disp(['det(R_sol) = ' num2str(det(R_sol))])
    disp(['resnorm = ' num2str(resnorm)])
    disp('R_sol = ') 
    disp(num2str(R_sol))
    %}
    %% 4.4 Subtract displacements due to the pure (theoretical) rotation
    xip_theo=(R_sol*xi0')';
    disp_theo=-xip_theo+xi0;       %= orig coords - rotn coords

    %% 4.5 Recreate original and cleaned displacement fields
    % Ux,Uy and Uz are displacements. Parameters without a prefix are
    % original and rb means rotated back. big_new_data_natural is the final
    % matrix similar to data corrected for rotation referenced to input x-y
    % frame.

    % rUi = disp_theo(:,1:3)  {theoretical displacments} 
    % Ui = cdata(:,4:6)       {measured displacements}
    % rbUi = Ui - rUi         {deformation displacements}

    rbUi = cdata(:,4:6) - disp_theo(:,1:3); 
    new_data=[cdata(:,1:3) rbUi];      %[CleanedCoords DeformDisps] N.B. cdata has been recentred c.f. tdata


else
    %if not RBD correction necessary
    new_data = cdata;
    eulerAngles = nan(1,3);
    rotCentre = nan(1,3);
    disp ' '
    disp(['RBD /mm = ' num2str(RBD,6)])
    disp('No correction for rotation requested');
end

%% 5. Recreate original and cleaned displacement fields
big_new_data_natural=data3D_col;

j=1;
for i=1:size(big_new_data_natural)
    if ~isnan(big_new_data_natural(i,6))
        big_new_data_natural(i,:)=new_data(j,:);
        j=j+1;    
    end
end
%% 7. Define outputs
eulerAngles = euler;
dataOut = big_new_data_natural;
disp(' ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Obsolete Code
%{
%% OC1   3.4.4 Extracting Euler angles from the rotation matrix
%{
t=zeros(3,1);

t(1)=atan2(R(2,3),R(3,3));
cos_t2=sqrt(R(1,1)^2+R(1,2)^2);
t(2)=atan2(-R(1,3),cos_t2);
t(3)=atan2(R(1,2),R(1,1));

disp('##ISSUE: check Euler angles and matrix generation##')

disp(['Euler angles /degrees: (alpha beta gamma) = (' num2str(radtodeg(t'),'% 1.3f ') ').'])

%% 3.4.5 Subtract displacements due to the pure (theoretical) rotation
c=cos(t);
s=sin(t);

% theoretical rotation matrix based on the calculated angles

R_theo=[c(2)*c(3), c(2)*s(3), -1*s(2); ...
     s(1)*s(2)*c(3)-c(1)*s(3), s(1)*s(2)*s(3)+c(1)*c(3), s(1)*c(2); ...
     c(1)*s(2)*c(3)+s(1)*s(3), c(1)*s(2)*s(3)-s(1)*c(3),c(1)*c(2)];

xip_theo=(R_theo*xi0')';
disp_theo=-xip_theo+xi0;       %= orig coords - rotn coords
%}
%%%%%%%%

%OC2
%{
%Redundant
%OC1
c=cos(euler);
s=sin(euler);

R_orig=[c(2)*c(3), s(1)*s(2)*c(3)-c(1)*s(3), c(1)*s(2)*c(3)+s(1)*s(3);...
        c(2)*s(3), s(1)*s(2)*s(3)+c(1)*c(3), c(1)*s(2)*s(3)-s(1)*c(3);...
        -s(2),     s(1)*c(2),                c(1)*c(2)              ];
%}
%%%%%%%%


%}
