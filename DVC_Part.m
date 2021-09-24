function outputname = DVC_Part(inO)
% the DVC whole directory data in .dat format
dvcdatapth = [inO '\3D Data Uxy.dat'];
% The .inp file where you got the nodes infromation
modelinp = [inO '\CrackedModel.inp'];%abaqus write input
% input file which will contain the BC from DVC DATA
outputname = [inO '\CrackedModel_BC.inp'];
vxsizDVC = 1; %Physical size of DVC = Pixel physical size* voxel size*(1-overlap)

%% Registration
[fedata,~] = readNodes(modelinp); Nodes = fedata(:,2:end);
OriPts = [max(Nodes(:,1)),max(Nodes(:,2)),min(Nodes(:,3));...
          max(Nodes(:,1)),max(Nodes(:,2)),max(Nodes(:,3)); ...
          min(Nodes(:,1)),min(Nodes(:,2)),min(Nodes(:,3))];
dvcdata = importdata(dvcdatapth);
dvcdata = dvcdata.data(:,1:6);   
% dvcdata(:,[2 5]) = -dvcdata(:,[2 5]);
DefPts = [max(dvcdata(:,1)),max(dvcdata(:,2)),min(dvcdata(:,3));...
          max(dvcdata(:,1)),max(dvcdata(:,2)),max(dvcdata(:,3)); ...
          min(dvcdata(:,1)),min(dvcdata(:,2)),min(dvcdata(:,3))];   
% transform dvc data to FE coordinates
dvcdata = AdjustmentBr0_v1beta(OriPts,DefPts,dvcdata);
% find forbidden nodes
radius = 10*vxsizDVC;
dvcdata = dvcdata(~(dvcdata(:,4)==0 & dvcdata(:,5)==0 & dvcdata(:,6)==0),:);
% interpolate DVC results on the FE mesh
fprintf('Training interpolators\n');
Fx = scatteredInterpolant(dvcdata(:,1),dvcdata(:,2),dvcdata(:,3),dvcdata(:,4),'natural');
Fy = scatteredInterpolant(dvcdata(:,1),dvcdata(:,2),dvcdata(:,3),dvcdata(:,5),'natural');
Fz = scatteredInterpolant(dvcdata(:,1),dvcdata(:,2),dvcdata(:,3),dvcdata(:,6),'natural');
fprintf('Interpolating Ux,Uy,Uz values\n');
fedata(:,5) = Fx(fedata(:,2),fedata(:,3),fedata(:,4));
fedata(:,6) = Fy(fedata(:,2),fedata(:,3),fedata(:,4));
fedata(:,7) = Fz(fedata(:,2),fedata(:,3),fedata(:,4));
fedata = fedata(sum(fedata(:,5:7),2)~=0,:);
%{
%% visualize the BC
inc = ceil(sqrt(length(dvcdata))/10); %adjust the number of arrows to be plotted
sc = 20; %adjust the length of the arrows
npts = size(fedata(1:inc:end,2),1);
figure, plot3(fedata(1:inc:end,2),fedata(1:inc:end,3),fedata(1:inc:end,4),...
    'kx','MarkerSize',6,'LineWidth',2);
xlabel('x [mm]');ylabel('y [mm]');zlabel('z [mm]') %x,y,z axes
hold on;h = quiver3(fedata(1:inc:end,2),fedata(1:inc:end,3),fedata(1:inc:end,4),...
    fedata(1:inc:end,5),fedata(1:inc:end,6),fedata(1:inc:end,7),sc);
set(h,'Color','r','LineWidth',2) %adjust the color, line width of the arrows
set(gca,'FontSize',18,'LineWidth',2,'box','on') %adjust the features of the axes
axis equal

figure, plot3(fedata(1:inc:end,2),fedata(1:inc:end,3),fedata(1:inc:end,4),'x');
xlabel('x');ylabel('y');zlabel('z')
hold on;quiver3(fedata(1:inc:end,2),fedata(1:inc:end,3),fedata(1:inc:end,4),...
    fedata(1:inc:end,5),fedata(1:inc:end,6),fedata(1:inc:end,7),sc)
%}
%% write the INPfile patches
fprintf('Writing INP patch\n');
% define the node sets (commented by Yixuan)
patchASSEM = cell(2*length(fedata),1);
for k=1:size(fedata,1)
    str1 = strcat('*Nset, nset=_SelimBC-',num2str(fedata(k,1)),...
        ', internal, instance=Part-1-1');
    str2 = strcat(num2str(fedata(k,1)),',');
    patchASSEM(k*2-1) =  cellstr(str1);
    patchASSEM(k*2) =  cellstr(str2);
%     fprintf(fileID2,'%s\n',str1);
%     fprintf(fileID2,'%s\n',str2);
end



patchBC = cell(length(fedata)+1,1);
patchBC(1) =  cellstr('*BOUNDARY, OP=NEW');
count=0;
for k=1:size(fedata,1)
    if ~isnan(fedata(k,5))
        str1 = strcat('_SelimBC-',num2str(fedata(k,1)),', 1, 1, ',num2str(fedata(k,5)));
        str2 = strcat('_SelimBC-',num2str(fedata(k,1)),', 2, 2, ',num2str(fedata(k,6)));
        str3 = strcat('_SelimBC-',num2str(fedata(k,1)),', 3, 3, ',num2str(fedata(k,7)));
        count = count+1;
        patchBC(count*3-1) =  cellstr(str1);
        patchBC(count*3) =  cellstr(str2);
        patchBC(count*3+1) =  cellstr(str3);
    end
end
% Parse blank INP file
fid = fopen(modelinp, 'r');
C = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

% Find where to inject the BC patch
injptBC = strfind(C{1}, '*Static');
rows = find(~cellfun('isempty', injptBC));
injptBC = rows +1;

%Find where to inject the ASSEM patch
injptASSEM = strfind(C{1}, '*End Assembly');
rows = find(~cellfun('isempty', injptASSEM));
injptASSEM = rows -1;

% Append data to the INP file
part1 = C{1}(1:injptASSEM);
part2 = C{1}(injptASSEM+1:injptBC);
part3 = C{1}(injptBC+1:end);
finalform = [part1;patchASSEM;part2;patchBC;part3];

% Write outputfile
fileID = fopen(outputname,'w');
for i=1:size(finalform,1)
    stri = finalform(i);
    fprintf(fileID,'%s\n',char(stri));
end
fclose(fileID);close all