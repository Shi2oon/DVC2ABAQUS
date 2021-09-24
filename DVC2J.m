function [offset,RadEulerAng,rotCentre,Abaqus,Len] = DVC2J(inO,Dir)
if strcmpi(Dir.Operation, 'DIC')
    [~,A] = fileparts(Dir.results);
    if contains(A,'3D_Integrated_Uxy')
        datum = load([erase(Dir.results,'.mat')  '.mat']);
        datum = [datum.M4.X(:)  datum.M4.Y(:) datum.M4.Z(:) datum.M4.Ux(:)...
                 datum.M4.Uy(:) datum.M4.Uz(:)];
    elseif exist([erase(Dir.results, '.mat') '.mat'],'file')
        datum = load([erase(Dir.results,'.mat')  '.mat']);     % save the file as all data
        datum = datum.alldata;
    elseif exist([erase(Dir.results, '.dat') '.dat'],'file')
        data  = importdata([erase(Dir.results, '.dat') '.dat']);
        datum = data.data;
    elseif exist(Dir.results,'file')
        data  = importdata(Dir.results);
        datum = data.data;
    end
%     if size(datum,2)==7
%         datum(datum(:,7)==0,1:6)=NaN;
%         datum = datum(:,1:6);
%     end
    datum = datum.*Dir.pixel_size;
    inO = fileparts(inO);
else 
%     [~,Dir.M4] = reshapeData(datum);
%     Dir.M4.X = Dir.M4.X1;   Dir.M4.Y = Dir.M4.Y1;   Dir.M4.Z = Dir.M4.Z1;else
    datum = [Dir.M4.X(:) 	Dir.M4.Y(:)     Dir.M4.Z(:)...
             Dir.M4.Ux(:)	Dir.M4.Uy(:)	Dir.M4.Uz(:)];
end
close all

%% CRACK COORDINATE (left crack)
if ~isfield(Dir,'xo');  Dir.xo = [];            Dir.yo = []; end
datum(isnan(datum))=0;
[Dir.xo,Dir.yo,Dir.M4] = get3DTip(datum,Dir.xo, Dir.yo);

%% get the paramters
Len = [min(Dir.M4.X(:))     max(Dir.M4.X(:));   % X dim	
       min(Dir.M4.Y(:))     max(Dir.M4.Y(:));   % Y dim
       max(Dir.M4.Z(:))     Dir.M4.MapStepSize;          % Z, Y crack coordinates
       Dir.xo(1)            Dir.xo(2);          % x crack coordinates
       Dir.yo(1)            Dir.yo(2);          % y crack coordinates
       Dir.M4.EdgeStepSize      ceil(length(unique(Dir.M4.X(:)))*0.5-2)]; % left hand
Len = round(Len,3);
Dir.M4.X = abs(Dir.M4.X-max(Dir.M4.X(:)));
[~,~, offset] = unist4Abaqus([],Dir.input_unit);
inO = [fileparts(inO) '\' Dir.unique];               mkdir(inO);
saveas(gcf, [inO '\DVC2ABAQUS Coodrinate.png']); close

%% create DVC .dat folder
tmp = sortrows([Dir.M4.X(:)     Dir.M4.Y(:)     Dir.M4.Z(:)  ...
                Dir.M4.Ux(:)    Dir.M4.Uy(:)    Dir.M4.Uz(:)],[3,1,2]);
[tmp,~ ] = reshapeData(tmp);
tmp(isnan(tmp))=0;
[~,dataum ] = reshapeData(tmp);
plotAllDis(dataum.Y1,dataum.X1,dataum.Z1,dataum.Ux,dataum.Uy,dataum.Uz,Dir.input_unit)
saveas(gcf,[inO '\Raw_DVC.fig']);    saveas(gcf,[inO '\Raw_DVC.tif']); close
        
inpFile = fopen([inO '\3D Data Uxy.dat'],'wt');
fprintf(inpFile, 'TITLE = "3D Data Uxy"\n');
fprintf(inpFile, 'VARIABLES = "x", "y", "z", "Vx", "Vy", "Vz", "isValid"\n');
l = [length(unique(Dir.M4.X(:))),length(unique(Dir.M4.Y(:))),...
    length(unique(Dir.M4.Z(:)))];
fprintf(inpFile, 'ZONE T="Frame 0", I=%d, J=%d, K=%d, F=POINT\n',l(1),l(2),l(3));
fprintf(inpFile, '%f %f %f %f %f %f %f\n',[Dir.M4.X(:)';Dir.M4.Y(:)';...
        Dir.M4.Z(:)'-min(unique(Dir.M4.Z(:)));Dir.M4.Ux(:)';Dir.M4.Uy(:)';...
        Dir.M4.Uz(:)';~isnan(Dir.M4.X(:))']);
fclose(inpFile);
%{
try
% remove rotation
[RadEulerAng,rotCentre,tmp] = shoemake_3D_v04_07_02_Abdo...
        (inO,'3D Data Uxy.dat',Dir.input_unit);
[~,dataum ] = reshapeData(tmp);
plotAllDis(dataum.Y1,dataum.X1,dataum.Z1,dataum.Ux,dataum.Uy,dataum.Uz,Dir.input_unit)
saveas(gcf,[inO '\Shoed_DVC.fig']);    saveas(gcf,[inO '\Shoed_DVC.tif']); close
l = [length(unique(dataum.X1(:))),length(unique(dataum.Y1(:))),...
    length(unique(dataum.Z1(:)))];
inpFile = fopen([inO '\3D Data Uxy.dat'],'wt');
fprintf(inpFile, 'ZONE T="Frame 0", I=%d, J=%d, K=%d, F=POINT\n',l(1),l(2),l(3));
fprintf(inpFile, '%f %f %f %f %f %f %f\n',[dataum.X1(:)';dataum.Y1(:)';...
        dataum.Z1(:)'-min(unique(dataum.Z1(:)));dataum.Ux(:)';dataum.Uy(:)';...
        dataum.Uz(:)';~isnan(dataum.X1(:))']);
fclose(inpFile);
catch err
    warning(['Rotation Removal failed because ' err.message])
    RadEulerAng =[];    rotCentre = [];
end
%}
RadEulerAng =[];        rotCentre = [];
%%
Greate3DCracked_v2(inO,round(Len,4),Dir,offset);  % Create a 3D crack in abaqus 
Abaqus = DVC_Part(inO);      % Merge the gotery with the boundray condtions
runBCJob(inO);      % Run the Job
printDVC(erase(Abaqus,'.inp'))
