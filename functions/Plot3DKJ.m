function [J,KI,KII,KIII,Der_Deg]=Plot3DKJ(saveto,offset,thickness,op)
set(0,'defaultAxesFontSize',25);       set(0,'DefaultLineMarkerSize',14)
if ~exist('thickness','var'); thickness = 1; end
if isempty(thickness);      thickness = 1; end
switch offset
    case 1e-6;  unit = '\mum';
    case 1e-3;  unit = 'mm';
    case 1;     unit = 'm';
    otherwise;  unit = 'unit';
end

%% read data
fid = fopen([erase(saveto, '.dat') '.dat'],'rt') ;
S = textscan(fid,'%s','Delimiter','\n');
S = S{1} ;      fclose('all');

%% DJ_DOWN
try
    idxS = strfind(S, 'CRACK-1_CRACK-1');
    idx1 = find(not(cellfun('isempty', idxS)));
    idxS = strfind(S, 'LABELS REFERENCED IN THE ABOVE TABLE');
    idx2 = find(not(cellfun('isempty', idxS)));
    J_in = S(idx1(3)+1:idx2(1)-3) ;% pick  nodes
catch
    idxS = strfind(S, 'H-OUTPUT-2_CRACK-1');
    idx1 = find(not(cellfun('isempty', idxS)));
    idxS = strfind(S, 'LABELS REFERENCED IN THE ABOVE TABLE');
    idx2 = find(not(cellfun('isempty', idxS)));
    J_in = S(idx1(3)+1:idx2(1)-3) ;
end

for iv = length(J_in):-1:1
    if isempty(J_in{iv,1})
        J_in(iv)=[];
    end
end

J_in = convertCharsToStrings(J_in);
countX=0; j=0;
for iv=1:size(J_in,1)
    %     Jvalue{iv,:} = cell2mat(cellfun(@str2num,J_in(iv),'UniformOutput',false));
    if isempty(J_in{iv})
        break;
    end
    Jvalue = textscan(J_in{iv},'%s');
    for ij = 1:length(Jvalue{1})
        if isnan(str2double(Jvalue{1}(ij)))
            j = j+1;countX=0;
        else
            countX = countX+1;
            J.Raw(j,countX) = str2double(Jvalue{1}(ij));
        end
    end
    clear Jvalue
end
% J.Raw = J.Raw(:,2:end);


%% KII and KII
try
    idxS = strfind(S, 'CRACK-2_CRACK-1');
    idx1 = find(not(cellfun('isempty', idxS)));
    K_in = S(idx1(3)+1:idx2(2)-3) ;
catch
    idxS = strfind(S, 'H-OUTPUT-3_CRACK-1');
    idx1 = find(not(cellfun('isempty', idxS)));
    K_in = S(idx1(3)+1:idx2(2)-3) ;
end

for iv = length(K_in):-1:1
    if isempty(K_in{iv,1})
        K_in(iv)=[];
    end
end
% driection
idxS = strfind(K_in, 'MTS   DIRECTION (DEG)');
idx1 = find(not(cellfun('isempty', idxS)));
if isempty(idx1)
    idxS = strfind(K_in, 'MERR  DIRECTION (DEG)');
    idx1 = find(not(cellfun('isempty', idxS)));
end
% KI,II,III
% K_in(idx1)=[];
B  = convertCharsToStrings(K_in);
countX=0;
for iv = 1:length(B)
    C{iv}  = textscan(B(iv),'%s');
    C{iv}  = C{iv}{1,1};
    format = 0;
    countX = countX+1;
    for ic = 1:length(C{iv})
        Conv = str2double(C{iv}(ic));
        if ~isnan(Conv)
            format = format+1;
            OutKJ(countX,format) = str2double(C{iv}(ic));
        end
    end
end
countX=1;
if ~isempty(idx1) 
    Step =5;
else
    Step =4;
end
for iv = 1:Step:size(OutKJ,1)
    KIRaw(countX:countX+size(OutKJ,2)-1)   = OutKJ(iv,:);
    KIIRaw(countX:countX+size(OutKJ,2)-1)  = OutKJ(iv+1,:);
    KIIIRaw(countX:countX+size(OutKJ,2)-1) = OutKJ(iv+2,:);
    if ~isempty(idx1)
        Der_Deg(countX:countX+size(OutKJ,2)-1) = OutKJ(iv+3,:);
        JKRaw(countX:countX+size(OutKJ,2)-1)  = OutKJ(iv+4,:);
    else
        Der_Deg = [];
        JKRaw(countX:countX+size(OutKJ,2)-1)  = OutKJ(iv+3,:);
    end
    countX = countX+size(OutKJ,2);
end

countX=1;
for iv = 1:size(J.Raw,1)
    for ij=1:size(J.Raw,2)
        KI.Raw(iv,ij)   = KIRaw(countX);
        KII.Raw(iv,ij)  = KIIRaw(countX);
        KIII.Raw(iv,ij) = KIIIRaw(countX);
        J.K.Raw(iv,ij)  = JKRaw(countX);
        countX=countX+1;
    end
end

%{
%% acessive countors
KI(KI==0)  =[];            KII(KII==0)=[];                     JK(JK==0)  =[];
KLM = OutJ(end);    [A,~] = ismember(OutJ,KLM);	  OutJ(A)= [];   OutJ(end+1) = KLM;
KLM = KI(end);      [A,~] = ismember(KI,KLM);	  KI(A)  = [];   KI(end+1)   = KLM;
KLM = KII(end);     [A,~] = ismember(KII,KLM);	  KII(A) = [];   KII(end+1)  = KLM;
KLM = JK(end);      [A,~] = ismember(JK,KLM);	  JK(A)  = [];   JK(end+1)   = KLM;

Kon = min([length(OutJ) length(JK) length(KI) length(KII)]);
OutJ(Kon:end)=[];   KI(Kon:end)=[]; KII(Kon:end)=[];    JK(Kon:end)=[];
%}
% J.K.Raw = J.K.Raw';     J.Raw = J.Raw';
% KI.Raw = KI.Raw';       KII.Raw = KII.Raw';
% KIII.Raw = KIII.Raw';
[a,b] = size(J.K.Raw);  A = ceil(a/2); B = ceil(b/2);
close all; plot(J.K.Raw(A,:)); text(1:b,J.K.Raw(A,:),string([1:b]))
set(gcf,'position',[98 311 1481 667])
oh = 11;%input('where to cut the contour? ');               close;
% oh=10; close
J.Raw = J.Raw(:,1:oh);      KI.Raw = KI.Raw(:,1:oh);  KII.Raw = KII.Raw(:,1:oh);
J.K.Raw = J.K.Raw(:,1:oh);  KIII.Raw = KIII.Raw(:,1:oh);        
[saveto,Ond] = fileparts(saveto);

%% plotting
J.Raw     = J.Raw./offset;             % in J/m^2
J.K.Raw   = J.K.Raw./offset;           % in J/m^2
KI.Raw    = KI.Raw./offset^1.5*1e-6;     % in MPa
KII.Raw   = KII.Raw./offset^1.5*1e-6;    % in MPa
KIII.Raw  = KIII.Raw./offset^1.5*1e-6;    % in MPa

contrs   = length(KI.Raw(1,:));        contrs = contrs - round(contrs*0.4);

% remove outliers
dic = real(ceil(-log10(nanmean(rmoutliers(J.K.Raw(contrs:end))))))+2;
if dic<1;   dic = 1; end
for iV=1:size(J.K.Raw,1)
    J.Ttrue(iV)   = round(mean(rmoutliers(J.Raw(iV,contrs:end))),dic);
    J.Tdiv(iV)    = round(std(rmoutliers(J.Raw(iV,contrs:end)),1),dic);
    J.K.Ttrue(iV) = round(mean(rmoutliers(J.K.Raw(iV,contrs:end))),dic);
    J.K.Tdiv(iV)  = round(std(rmoutliers(J.K.Raw(iV,contrs:end)),1),dic);
    KIII.Ttrue(iV)= round(mean(rmoutliers(KIII.Raw(iV,contrs:end))),dic);
    KIII.Tdiv(iV) = round(std(rmoutliers(KIII.Raw(iV,contrs:end)),1),dic);
    KI.Ttrue(iV)  = round(mean(rmoutliers(KI.Raw(iV,contrs:end))),dic);
    KI.Tdiv(iV)   = round(std(rmoutliers(KI.Raw(iV,contrs:end)),1),dic);
    KII.Ttrue(iV) = round(mean(rmoutliers(KII.Raw(iV,contrs:end))),dic);
    KII.Tdiv(iV)  = round(std(rmoutliers(KII.Raw(iV,contrs:end)),1),dic);
end

%%
thickness = [1:iV]*thickness;
    J.true   = round(mean(rmoutliers(mean(rmoutliers(J.Raw(:,contrs:end))))),dic);
    J.div    = round(std(rmoutliers(mean(rmoutliers(J.Raw(:,contrs:end)))),1),dic);
    J.K.true = round(mean(rmoutliers(mean(rmoutliers(J.K.Raw(:,contrs:end))))),dic);
    J.K.div  = round(std(rmoutliers(mean(rmoutliers(J.K.Raw(:,contrs:end)))),1),dic);
    KIII.true  = round(mean(rmoutliers(mean(rmoutliers(KIII.Raw(:,contrs:end))))),dic);
    KIII.div   = round(std(rmoutliers(mean(rmoutliers(KIII.Raw(:,contrs:end)))),1),dic);
    KI.true  = round(mean(rmoutliers(mean(rmoutliers(KI.Raw(:,contrs:end))))),dic);
    KI.div   = round(std(rmoutliers(mean(rmoutliers(KI.Raw(:,contrs:end)))),1),dic);
    KII.true = round(mean(rmoutliers(mean(rmoutliers(KII.Raw(:,contrs:end))))),dic);
    KII.div  = round(std(rmoutliers(mean(rmoutliers(KII.Raw(:,contrs:end)))),1),dic);

if size(J.K.Raw,1)<4
    for iV=1:size(J.K.Raw,1)
        fig=figure;set(fig,'defaultAxesColorOrder',[[0.4 0 0.5]; [0 0 0]]);
        yyaxis left;
        plot(KI.Raw(iV,:),'b--d','MarkerEdgeColor','b','LineWidth',2); hold on;
        plot(KII.Raw(iV,:),'r--s','MarkerEdgeColor','r','LineWidth',2);
        plot(KIII.Raw(iV,:),'g--<','MarkerEdgeColor','g','LineWidth',2); hold off
        ylabel('K (MPa m^{0.5})');
%         if min([KII.true(:); KI.true(:); KIII.true(:)])>0;     ylim([0 inf]);      end
        yyaxis right;
        plot(J.K.Raw(iV,:),'k--o','MarkerEdgeColor','k','LineWidth',2,'MarkerFaceColor','k');
        ylabel('J [J/m^2]');axis tight
        xlabel('Contour Number');
        legend(['K_{I} = ' num2str(KI.Ttrue(iV)) ' ± ' num2str(KI.Tdiv(iV) ) ' MPa\surdm' ],...
            ['K_{II} = ' num2str(KII.Ttrue(iV)) ' ± ' num2str(KII.Tdiv(iV) ) ' MPa\surdm' ],...
            ['K_{III} = ' num2str(KIII.Ttrue(iV)) ' ± ' num2str(KIII.Tdiv(iV) ) ' MPa\surdm' ],...
            ['J_{integral} = ' num2str(J.K.Ttrue(iV)) ' ± ' ...
            num2str(J.K.Tdiv(iV) ) ' J/m^2'],'location','northoutside','box','off');
        set(gcf,'WindowStyle','normal');
        set(gcf,'position',[60,10,850,990]);
        axis tight
        box off; saveas(gcf, [saveto '\' Ond '_' num2str(iV) '_KI, KII and J.fig']);
        saveas(gcf, [saveto '\' Ond '_' num2str(iV) '_KI, KII and J.tif']);    close
    end
end
if exist('op','var')
    for iV=op
        fig=figure;set(fig,'defaultAxesColorOrder',[[0.4 0 0.5]; [0 0 0]]);
        yyaxis left;
        plot(KI.Raw(iV,:),'b--d','MarkerEdgeColor','b','LineWidth',2); hold on;
        plot(KII.Raw(iV,:),'r--s','MarkerEdgeColor','r','LineWidth',2);
        plot(KIII.Raw(iV,:),'g--<','MarkerEdgeColor','g','LineWidth',2); hold off
        ylabel('K (MPa m^{0.5})');
%         if min([KII.true(:); KI.true(:); KIII.true(:)])>0;     ylim([0 inf]);      end
        yyaxis right;
        plot(J.K.Raw(iV,:),'k--o','MarkerEdgeColor','k','LineWidth',2,'MarkerFaceColor','k');
        ylabel('J [J/m^2]');axis tight
        xlabel('Contour Number');
        legend(['K_{I} = ' num2str(KI.Ttrue(iV)) ' ± ' num2str(KI.Tdiv(iV) ) ' MPa\surdm' ],...
            ['K_{II} = ' num2str(KII.Ttrue(iV)) ' ± ' num2str(KII.Tdiv(iV) ) ' MPa\surdm' ],...
            ['K_{III} = ' num2str(KIII.Ttrue(iV)) ' ± ' num2str(KIII.Tdiv(iV) ) ' MPa\surdm' ],...
            ['J_{integral} = ' num2str(J.K.Ttrue(iV)) ' ± ' ...
            num2str(J.K.Tdiv(iV) ) ' J/m^2'],'location','northoutside','box','off');
        set(gcf,'WindowStyle','normal');
        set(gcf,'position',[60,10,850,990]);
        axis tight
        box off; saveas(gcf, [saveto '\' Ond '_' num2str(iV) '_KI, KII and J.fig']);
        saveas(gcf, [saveto '\' Ond '_' num2str(iV) '_KI, KII and J.tif']);    close
    end
end

%% VALUES PER THICKNESS Raw
C = unique(diff(thickness(:)));   
fig=figure;set(fig,'defaultAxesColorOrder',[[0.4 0 0.5]; [0 0 0]]); yyaxis left;
errorbar(thickness,KI.Ttrue,KI.Tdiv,'b--d','MarkerEdgeColor','b',...
    'LineWidth',2,'markersize',10); hold on;
errorbar(thickness,KII.Ttrue,KII.Tdiv,'r--s','MarkerEdgeColor','r',...
    'LineWidth',2,'markersize',10);
errorbar(thickness,KIII.Ttrue,KIII.Tdiv,'g--<','MarkerEdgeColor','g',...
    'LineWidth',2,'markersize',10); hold off
if min([KII.Ttrue(:); KI.Ttrue(:); KIII.Ttrue(:)])>0;     ylim([0 inf]);      end
ylabel('K (MPa m^{0.5})');      yyaxis right;
errorbar(thickness,J.K.Ttrue,J.K.Tdiv,'k--o','MarkerEdgeColor','k','LineWidth',2,...
    'MarkerFaceColor','k','markersize',10);
ylabel('J [J/m^2]');        ylim([0 inf]);
xlabel(['Thickness (' unit ')']);
legend(['K_{I} = ' num2str(KI.true) ' ± ' num2str(KI.div) ' MPa\surdm' ],...
    ['K_{II} = ' num2str(KII.true) ' ± ' num2str(KII.div) ' MPa\surdm' ],...
    ['K_{III} = ' num2str(KIII.true) ' ± ' num2str(KIII.div) ' MPa\surdm' ],...
    ['J_{integral} = ' num2str(J.K.true) ' ± ' ...
    num2str(J.K.div) ' J/m^2'],'location','northoutside','box','off');
set(gcf,'WindowStyle','normal'); set(gcf,'position',[60,10,850,990]);
xlim([0 max(thickness)+C(1)])
box off; saveas(gcf, [saveto '\' Ond '_Thickness_KI, KII and J_Raw.fig']);
saveas(gcf, [saveto '\' Ond '_Thickness_KI, KII and J_Raw.tif']);    close

% trimmed
[~,B]=rmoutliers(J.K.Ttrue);      thickness(B)=NaN;
fig=figure;set(fig,'defaultAxesColorOrder',[[0.4 0 0.5]; [0 0 0]]); yyaxis left;
errorbar(thickness,KI.Ttrue,KI.Tdiv,'b--d','MarkerEdgeColor','b',...
    'LineWidth',2,'markersize',10); hold on;
errorbar(thickness,KII.Ttrue,KII.Tdiv,'r--s','MarkerEdgeColor','r',...
    'LineWidth',2,'markersize',10);
errorbar(thickness,KIII.Ttrue,KIII.Tdiv,'g--<','MarkerEdgeColor','g',...
    'LineWidth',2,'markersize',10); hold off
if min([KII.Ttrue(:); KI.Ttrue(:); KIII.Ttrue(:)])>0;     ylim([0 inf]);      end
ylabel('K (MPa m^{0.5})');      yyaxis right;
errorbar(thickness,J.K.Ttrue,J.K.Tdiv,'k--o','MarkerEdgeColor','k','LineWidth',2,...
    'MarkerFaceColor','k','markersize',10);
ylabel('J [J/m^2]');        ylim([0 inf]);
xlabel(['Thickness (' unit ')']);
legend(['K_{I} = ' num2str(KI.true) ' ± ' num2str(KI.div) ' MPa\surdm' ],...
    ['K_{II} = ' num2str(KII.true) ' ± ' num2str(KII.div) ' MPa\surdm' ],...
    ['K_{III} = ' num2str(KIII.true) ' ± ' num2str(KIII.div) ' MPa\surdm' ],...
    ['J_{integral} = ' num2str(J.K.true) ' ± ' ...
    num2str(J.K.div) ' J/m^2'],'location','northoutside','box','off');
set(gcf,'WindowStyle','normal'); set(gcf,'position',[60,10,850,990]);
xlim([0 max(thickness)+C(1)])
box off; saveas(gcf, [saveto '\' Ond '_Thickness_KI, KII and J.fig']);
saveas(gcf, [saveto '\' Ond '_Thickness_KI, KII and J.tif']);    close

%{
subplot(1,2,1)
errorbar(thickness,J.K.Ttrue,J.K.Tdiv,'k--o','MarkerEdgeColor','k',...
            'LineWidth',2,'MarkerFaceColor','k','markersize',10); hold on;
errorbar(thickness,J.K.Ttrue,J.K.Tdiv,'r--s','MarkerEdgeColor',...
            'r','LineWidth',2,'markersize',10); hold off
xlabel(['Thickness (' unit ')']);   xlim([0 max(thickness)+abs(thickness(1)-thickness(2))])
ylabel('J [J/m^2]');
legend('J_{integral}', 'J_{SIF}','location','northoutside','box','off');
subplot(1,2,2)
plot(J.Avg,'k--o','MarkerEdgeColor',...
            'k','LineWidth',2,'MarkerFaceColor','k','markersize',10); hold on;
plot(J.K.Avg,'r--s','MarkerEdgeColor',...
            'r','LineWidth',2,'markersize',10); hold off
xlabel('Contour Number');           xlim([0 length(KI.Avg)+1])
ylabel('J [J/m^2]');
legend(['J_{integral} = ' num2str(J.true) ' ± ' num2str(J.div) ' J/m^2'],...
       ['J_{SIF} = ' num2str(J.K.true) ' ± ' num2str(J.div) ' J/m^2'],...
       'location','northoutside','box','off');
set(gcf,'WindowStyle','normal');    set(gcf,'position',[60 90 1298 882]);
box off; saveas(gcf, [saveto '\' Ond '_Thickness_KI, KII and J.fig']);
saveas(gcf, [saveto '\' Ond '_Thickness_KI, KII and J.tif']);    close
%}
% save([saveto '\' Ond],'J','KI','KII','KIII')
%%
if ~exist(fullfile(saveto,[Ond '_S12N.tif']),'file')
    fileid = fullfile(saveto, 'Abaqus_script_3.py');
    fileID = fopen(fileid,'w');
    ImportModules(fileID,saveto)
    SaveReport(fileID,saveto,Ond)          %% save report
    printTiffs(fileID,saveto,'N',Ond,3);   %% save tiff images for U and S
    printTiffs(fileID,saveto,'D',Ond,3);   %% save tiff images for deformaed U and S
    fclose(fileID);
    PWD =pwd; cd(saveto);
    system(['abaqus cae ','noGUI','=Abaqus_script_3.py']); % Windows system
    cd(PWD);
end

