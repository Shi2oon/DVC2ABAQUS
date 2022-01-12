% create a 30 MPa mode-I field (for mode II input 'II' for mixed mode input 'fun'
% then decide number of elements LxL. here I choose L = 8, to get faster
% results
restoredefaultpath;clc;clear;close all
addpath(genpath([pwd '\functions']));
%                                           [KI,KII,KIII],'fun',L
[files,Operation,unit,stif] = Westergaard_3D([30,10,50],'fun',25);
MatP.Operation = Operation{1,2};
MatP.input_unit   = unit{1,2};        % meter (m) or milmeter (mm) or micrometer(um);
MatP.pixel_size   = 1;           % if DIC values are in pixel, 1 if in physical units;
MatP.results = files{1,2};
MatP.stressstat = Operation{1,2}; % plane stress or strain
MatP.Mat = 'ferrite'; % materil name (any will do)
MatP.type = 'E'; % E for istropic linear elastic, A for anistropic with
% stinfess tensor as input and R for elastic plastic mateirals
MatP.E = 210e9;  % young's modulus
MatP.nu = 0.3;   % possiov's ratio
MatP.unique = 'calibration'; % unique name

%% see the error with the crack location in x and y (z not applicable);
[offset,RadEulerAng,rotCentre,Abaqus,Len] = DVC2J(MatP);% input to abaqus
% Plot results
[J,KI,KII,KIII,Direction] = Plot3DKJ(erase(Abaqus,'.inp'),offset,Len(3,1));

%{
%%  more details INPUT MATERIAL PROPERTIES AND DATA
Dir.type  = 'M'; % for engineering constants'E'
Dir.E = [145.3, 7.115,     7.115]*1e9;
Dir.nu = [0.256, 0.256, 0.48];
Dir.G = [4.7, 4.7, 3.0]*1e9;


Dir.type  = 'R'; % Ramberg-Osgood
Dir.Exponent = 26.67;
Dir.Yield_offset = 1.24;
Dir.yield = 4E9;			%Yield Stress [Pa]


Dir.type  = 'A'; % Elastic-Anisotropic
Dir.Stiffness = []; % a 6*6 matrix
%}
