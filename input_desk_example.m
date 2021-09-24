clc;clear;close all
% create a 30 MPa mode-I field (for mode II input 'II' for mixed mode input 'fun'
% then decide number of elements LxL. here I choose L = 8, to get faster
% results
clc;clear;close all
[files,Operation,unit,stif] = Westergaard_3D(30,'I',20);
MatP.Operation = Operation{1,2};
MatP.input_unit   = unit{1,2};        % meter (m) or milmeter (mm) or micrometer(um);
MatP.pixel_size   = 1;           % if DIC values are in pixel, 1 if in physical units;
MatP.results = files{1,2};
MatP.stressstat = Operation{1,2}; % plane stress or strain
MatP.Mat = 'ferrite'; % materil name (any will do)
MatP.type = 'E'; % E for istropic linear elastic, A for anistropic with 
                 % stinfess tensor as input and R for elastic plastic mateirals
MatP.E = 220e9;  % young's modulus
MatP.nu = 0.3;   % possiov's ratio
MatP.unique = 'calibration'; % unique name

% input to abaqus
[offset,RadEulerAng,rotCentre,Abaqus,Len] = DVC2J(MatP);

% Plot results
[J,KI,KII,KIII,Direction] = Plot3DKJ(erase(Abaqus,'.inp'),offset,Len(3,1));