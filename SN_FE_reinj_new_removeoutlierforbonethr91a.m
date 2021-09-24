clear
clc

dvcdatapth = 'W:\bone\thr61bslices\bone16bitpreloadvsload1-DVC_Vec_64x64x64_75ov-20=unknownfraction-0.8\RBM corrected/text-All.txt';
fe_pth = 'W:\bone\thr61bslices/node.xlsx';
modelinp = 'W:\bone\thr61bslices/Job-1.inp';
outputname = 'W:\bone\thr61bslices\bone16bitpreloadvsload1-DVC_Vec_64x64x64_75ov-20=unknownfraction-0.8\RBM corrected/PatchedInp.inp';

% Registration
OriPts = [2.92813 0.259305 3.806792;2.96532 2.79537 3.806792;0.713457 3.427302 3.933552]; %FE
DefPts=[0.2266068 -2.524738 0;0.288603 -3.104082 3.933552;1.220684 -0.1902607 3.933552]; %DVC

% import dvc data
dvcdata = importdata(dvcdatapth);
dvcdata = dvcdata.data;
% put dvc to right handed coordi
dvcdata(:,[2 5]) = -dvcdata(:,[2 5]);

% transform dvc data to FE coordinates
dvcdata = AdjustmentBr0_v1beta(OriPts,DefPts,dvcdata);

% import FE nodes position
[fedata,~,~] = xlsread(fe_pth,'Sheet1');
fedata = fedata(:,1:4);

% find forbidden nodes
% radius = 32*12*3.25*1E-3;
% forbidnod = fedata(fedata(:,3)<3.4 & fedata(:,4)>3.422-radius & fedata(:,4) <3.422+radius,1);
% forbidnod = [forbidnod; fedata(fedata(:,3)>5.0)];
% forbidnod = unique(forbidnod);

% Remove forbidden nodes
% [~,authlines] = setdiff(fedata(:,1),forbidnod);
% fedata = fedata(authlines,:);
figure, plot3(fedata(1:10:end,2),fedata(1:10:end,3),fedata(1:10:end,4),'x');

% Clean DVC data
dvcdata = dvcdata(~(dvcdata(:,4)==0 & dvcdata(:,5)==0 & dvcdata(:,6)==0),:);
dvcdata(dvcdata(:,3)<0.06,:) = [];

% Try to remove outliers
%dvcdata(dvcdata(:,4)<-0.0202 | dvcdata(:,5)<-0.0182 | dvcdata(:,6)<-0.0158,:)=[];
dvcdata(dvcdata(:,4)>0.018 | dvcdata(:,5)>0.017 | dvcdata(:,6)>0.011,:)=[];

dvcdata(dvcdata(:,1)<min(fedata(:,2)),:) = []; dvcdata(dvcdata(:,2)<min(fedata(:,3)),:) = []; dvcdata(dvcdata(:,3)<min(fedata(:,4)),:) = [];
dvcdata(dvcdata(:,1)>max(fedata(:,2)),:) = []; dvcdata(dvcdata(:,2)>max(fedata(:,3)),:) = []; dvcdata(dvcdata(:,3)>max(fedata(:,4)),:) = [];
fedata(fedata(:,2)<min(dvcdata(:,1)),:) = []; fedata(fedata(:,3)<min(dvcdata(:,2)),:) = [];fedata(fedata(:,4)<min(dvcdata(:,3)),:) = [];
fedata(fedata(:,2)>max(dvcdata(:,1)),:) = []; fedata(fedata(:,3)>max(dvcdata(:,2)),:) = [];fedata(fedata(:,4)>max(dvcdata(:,3)),:) = [];

% interpolate DVC results on the FE mesh
fprintf('Training interpolators\n');
Fx = scatteredInterpolant(dvcdata(:,1),dvcdata(:,2),dvcdata(:,3),dvcdata(:,4),'natural');
Fy = scatteredInterpolant(dvcdata(:,1),dvcdata(:,2),dvcdata(:,3),dvcdata(:,5),'natural');
Fz = scatteredInterpolant(dvcdata(:,1),dvcdata(:,2),dvcdata(:,3),dvcdata(:,6),'natural');

fprintf('Interpolating Ux,Uy,Uz values\n');
fedata(:,5) = Fx(fedata(:,2),fedata(:,3),fedata(:,4));
fedata(:,6) = Fy(fedata(:,2),fedata(:,3),fedata(:,4));
fedata(:,7) = Fz(fedata(:,2),fedata(:,3),fedata(:,4));

%remove points where BC are set to 0,0,0
fedata = fedata(sum(fedata(:,5:7),2)~=0,:);

% write the INPfile patches
fprintf('Writing INP patch\n');

% fileID = fopen('./PATCHS/NEWpatchBC.txt','w');
% fprintf(fileID,'%s\n','*BOUNDARY, OP=NEW');
patchBC = cell(3*length(fedata)+1,1);
patchBC(1) =  cellstr('*BOUNDARY, OP=NEW');
for k=1:length(fedata)
    str1 = strcat('_SelimBC-',num2str(fedata(k,1)),', 1, 1, ',num2str(fedata(k,5)));
    str2 = strcat('_SelimBC-',num2str(fedata(k,1)),', 2, 2, ',num2str(fedata(k,6)));
    str3 = strcat('_SelimBC-',num2str(fedata(k,1)),', 3, 3, ',num2str(fedata(k,7)));
    patchBC(k*3-1) =  cellstr(str1);
    patchBC(k*3) =  cellstr(str2);
    patchBC(k*3+1) =  cellstr(str3);
%     fprintf(fileID,'%s\n',str1);
%     fprintf(fileID,'%s\n',str2);
%     fprintf(fileID,'%s\n',str3);
end
% fclose(fileID);

% fileID2 = fopen('./PATCHS/NEWpatchASEM.txt','w');
patchASSEM = cell(2*length(fedata),1);
for k=1:length(fedata)
    str1 = strcat('*Nset, nset=_SelimBC-',num2str(fedata(k,1)),', internal, instance=Part-1-1');
    str2 = strcat(num2str(fedata(k,1)),',');
    patchASSEM(k*2-1) =  cellstr(str1);
    patchASSEM(k*2) =  cellstr(str2);
%     fprintf(fileID2,'%s\n',str1);
%     fprintf(fileID2,'%s\n',str2);
end
% fclose(fileID2);

% Parse blank INP file
fid = fopen(modelinp, 'r');
C = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

%Find where to inject the BC patch
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
for i=1:length(finalform)
    stri = finalform(i);
    fprintf(fileID,'%s\n',char(stri));
end
fclose(fileID);





