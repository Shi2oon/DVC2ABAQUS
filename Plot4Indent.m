function Plot4Indent(saveto,offset,thickness,nY)
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
idxS = strfind(S, ['CRJ' num2str(nY) '_CRACK-' num2str(nY)]);
idx1 = find(not(cellfun('isempty', idxS)));
idxS = strfind(S, 'LABELS REFERENCED IN THE ABOVE TABLE');
idx2 = find(not(cellfun('isempty', idxS)));
J_in = S(idx1(3)+1:idx2(nY)-3) ;

J_in = convertCharsToStrings(J_in);
countX=0; j=1;
for iv=1:size(J_in,1)
    %     Jvalue{iv,:} = cell2mat(cellfun(@str2num,J_in(iv),'UniformOutput',false));
    if isempty(J_in{iv})
        j = j+1;   countX=0;
    else
        Jvalue = textscan(J_in{iv},'%s');
        for ij = 1:length(Jvalue{1})
            %         if isnan(str2double(Jvalue{1}(ij)))
            %             j = j+1;countX=0;
            %         else
            countX = countX+1;
            J.Raw(j,countX) = str2double(Jvalue{1}(ij));
            %         end
        end
        clear Jvalue
    end
end
J.Raw(J.Raw==0)=NaN;
J.Raw(:,1)=[];

%% KII and KII
idxS = strfind(S, ['CRK' num2str(nY) '_CRACK-' num2str(nY)]);
idx1 = find(not(cellfun('isempty', idxS)));
K_in = S(idx1(3)+1:idx2(4+nY)-3) ;

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
if ~isempty(idx1)
    DIRECTION = convertCharsToStrings(K_in(idx1));
    countX=0;
    for iv=1:size(DIRECTION,1)
        O{iv}  = textscan(DIRECTION(iv),'%s');
        O{iv}  = O{iv}{1,1};
        format = 0;
        countX = countX+1;
        for ic = 1:length(O{iv})
            Conv = str2double(O{iv}(ic));
            if isnan(Conv)==0
                format = format+1;
                Der_Deg(countX,format) = str2double(O{iv}(ic));
            end
        end
    end
    Der_Deg = Der_Deg';
    Der_Deg = Der_Deg(:);
else
    Der_Deg = [];
end

K_in(idx1)=[];
B  = convertCharsToStrings(K_in);
countX=0;
for iv = 1:length(B)
    C{iv}  = textscan(B(iv),'%s');
    C{iv}  = C{iv}{1,1};
    format = 0;
    countX = countX+1;
    for ic = 1:length(C{iv})
        Conv = str2double(C{iv}(ic));
        if isnan(Conv)==0
            format = format+1;
            OutKJ(countX,format) = str2double(C{iv}(ic));
        end
    end
end
countX=1;
for iv = 1:4:size(OutKJ,1)
    KIRaw(countX:countX+size(OutKJ,2)-1)   = OutKJ(iv,:);
    KIIRaw(countX:countX+size(OutKJ,2)-1)  = OutKJ(iv+1,:);
    KIIIRaw(countX:countX+size(OutKJ,2)-1) = OutKJ(iv+2,:);
    JKRaw(countX:countX+size(OutKJ,2)-1)  = abs(OutKJ(iv+3,:));
    countX = countX+size(OutKJ,2);
end

countX=1;
for iv = 1:size(J.Raw,1)
    for ij=1:size(J.Raw,2)
        if ~isnan(J.Raw(iv,ij))
            KI.Raw(iv,ij)   = KIRaw(countX);
            KII.Raw(iv,ij)  = KIIRaw(countX);
            KIII.Raw(iv,ij) = KIIIRaw(countX);
            J.K.Raw(iv,ij)  = JKRaw(countX);
            countX=countX+1;
        end
    end
%     countX=countX+1;
end

[a,b] = size(KI.Raw);  A = ceil(a/2); B = ceil(b/2);
close all; plot(J.Raw(A,:)); text(1:b,J.Raw(A,:),string([1:b]))
set(gcf,'position',[98 311 1481 667])
oh = input('where to cut the contour? ');               close;
J.Raw = J.Raw(:,1:oh);      KI.Raw = KI.Raw(:,1:oh);  KII.Raw = KII.Raw(:,1:oh);
J.K.Raw = J.K.Raw(:,1:oh);  KIII.Raw = KIII.Raw(:,1:oh);
if size(KI.Raw,1)>2
    plot(J.Raw(:,B)); text(1:a,J.Raw(:,B),string([1:a]))
    set(gcf,'position',[98 311 1481 667])
    
    oh = input('where to cut the thickness [beign end]? '); close;
    J.Raw = J.Raw(oh(1):oh(2),:);       KI.Raw = KI.Raw(oh(1):oh(2),:);
    KII.Raw = KII.Raw(oh(1):oh(2),:);   J.K.Raw = J.K.Raw(oh(1):oh(2),:);
    KIII.Raw = KIII.Raw(oh(1):oh(2),:);
else
    J.Raw = J.Raw(1:size(KI.Raw,1),:);
end;        [saveto,ok] = fileparts(saveto);
Ond = ['Cr_' num2str(nY)];

%% plotting
J.Raw     = J.Raw./offset;             % in J/m^2
J.K.Raw   = J.K.Raw./offset;           % in J/m^2
KI.Raw    = KI.Raw./offset^1.5*1e-6;     % in MPa
KII.Raw   = KII.Raw./offset^1.5*1e-6;    % in MPa
KIII.Raw  = KIII.Raw./offset^1.5*1e-6;    % in MPa

J.Avg     = mean(J.Raw,1);             % in J/m^2
J.K.Avg   = mean(J.K.Raw,1);          % in J/m^2
KI.Avg    = mean(KI.Raw,1);     % in MPa
KII.Avg   = mean(KII.Raw,1);    % in MPa
KIII.Avg  = mean(KIII.Raw,1);  % in MPa

contrs   = length(KI.Avg);        contrs = contrs - round(contrs*0.4);

% remove outliers
dic = real(ceil(-log10(nanmean(rmoutliers(J.Raw(contrs:end))))))+2;
if dic<1;   dic = 1; end
J.true   = round(mean(rmoutliers(J.Avg(contrs:end))),dic);
J.div    = round(std(rmoutliers(J.Avg(contrs:end)),1),dic);
J.K.true = round(mean(rmoutliers(J.K.Avg(contrs:end))),dic);
J.K.div  = round(std(rmoutliers(J.K.Avg(contrs:end)),1),dic);
KIII.true  = round(mean(rmoutliers(KIII.Avg(contrs:end))),dic);
KIII.div   = round(std(rmoutliers(KIII.Avg(contrs:end)),1),dic);
KI.true  = round(mean(rmoutliers(KI.Avg(contrs:end))),dic);
KI.div   = round(std(rmoutliers(KI.Avg(contrs:end)),1),dic);
KII.true = round(mean(rmoutliers(KII.Avg(contrs:end))),dic);
KII.div  = round(std(rmoutliers(KII.Avg(contrs:end)),1),dic);

fig=figure;set(fig,'defaultAxesColorOrder',[[0.4 0 0.5]; [0 0 0]]);
yyaxis left;
plot(KI.Avg,'b--d','MarkerEdgeColor','b','LineWidth',2); hold on;
plot(KII.Avg,'r--s','MarkerEdgeColor','r','LineWidth',2);
plot(KIII.Avg,'g--<','MarkerEdgeColor','g','LineWidth',2); hold off
ylabel('K (MPa m^{0.5})');
if min([KII.Avg(:); KI.Avg(:); KIII.Avg(:)])>0;     ylim([0 inf]);      end
yyaxis right;
plot(J.Avg,'k--o','MarkerEdgeColor','k','LineWidth',2,'MarkerFaceColor','k');
ylabel('J [J/m^2]');
xlabel('Contour Number');
legend(['K_{I} = ' num2str(KI.true) ' ± ' num2str(KI.div) ' MPa\surdm' ],...
    ['K_{II} = ' num2str(KII.true) ' ± ' num2str(KII.div) ' MPa\surdm' ],...
    ['K_{III} = ' num2str(KIII.true) ' ± ' num2str(KIII.div) ' MPa\surdm' ],...
    ['J_{integral} = ' num2str(J.true) ' ± ' ...
    num2str(J.div) ' J/m^2'],'location','northoutside','box','off');
set(gcf,'WindowStyle','normal');
set(gcf,'position',[60,10,850,990]); xlim([0 length(KI.Avg)+1])
box off; saveas(gcf, [saveto '\' Ond '_Avg_KI, KII and J.fig']);
saveas(gcf, [saveto '\' Ond '_Avg_KI, KII and J.tif']);    close

fig=figure;set(fig,'defaultAxesColorOrder',[[0.4 0 0.5]; [0 0 0]]);
yyaxis left;
plot(KI.Avg,'b--d','MarkerEdgeColor','b','LineWidth',2); hold on;
plot(KII.Avg,'r--s','MarkerEdgeColor','r','LineWidth',2);
plot(KIII.Avg,'g--<','MarkerEdgeColor','g','LineWidth',2); hold off
ylabel('K (MPa m^{0.5})');
if min([KII.Avg(:); KI.Avg(:); KIII.Avg(:)])>0;     ylim([0 inf]);      end
yyaxis right;
plot(J.K.Avg,'k--o','MarkerEdgeColor','k','LineWidth',2,'MarkerFaceColor','k');
ylabel('J [J/m^2]');
xlabel('Contour Number');
legend(['K_{I} = ' num2str(KI.true) ' ± ' num2str(KI.div) ' MPa\surdm' ],...
    ['K_{II} = ' num2str(KII.true) ' ± ' num2str(KII.div) ' MPa\surdm' ],...
    ['K_{III} = ' num2str(KIII.true) ' ± ' num2str(KIII.div) ' MPa\surdm' ],...
    ['J_{integral} = ' num2str(J.K.true) ' ± ' ...
    num2str(J.K.div) ' J/m^2'],'location','northoutside','box','off');
set(gcf,'WindowStyle','normal');
set(gcf,'position',[60,10,850,990]); xlim([0 length(KI.Avg)+1])
box off; saveas(gcf, [saveto '\' Ond '_Avg_KI, KII and JK.fig']);
saveas(gcf, [saveto '\' Ond '_Avg_KI, KII and JK.tif']);    close

plot(KI.Avg,'b--d','MarkerEdgeColor','b','LineWidth',2); hold on;
plot(KII.Avg,'r--s','MarkerEdgeColor','r','LineWidth',2);
plot(KIII.Avg,'g--<','MarkerEdgeColor','g','LineWidth',2); hold off
ylabel('K (MPa m^{0.5})');
if min([KII.Avg(:); KI.Avg(:); KIII.Avg(:)])>0;     ylim([0 inf]);      end
xlabel('Contour Number');
legend(['K_{I} = ' num2str(KI.true) ' ± ' num2str(KI.div) ' MPa\surdm' ],...
    ['K_{II} = ' num2str(KII.true) ' ± ' num2str(KII.div) ' MPa\surdm' ],...
    ['K_{III} = ' num2str(KIII.true) ' ± ' num2str(KIII.div) ' MPa\surdm' ],...
    'location','northoutside','box','off');
set(gcf,'WindowStyle','normal');
set(gcf,'position',[60,10,850,990]); xlim([0 length(KI.Avg)+1])
box off; saveas(gcf, [saveto '\' Ond '_Avg_KI, KII and JK.fig']);
saveas(gcf, [saveto '\' Ond '_Avg_KI-II-III.tif']);
saveas(gcf, [saveto '\' Ond '_Avg_KI-II-III.fig']); close
%%
for iV=1:size(KI.Raw,1)
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
thickness = [1:iV]*thickness;

for iV=1:size(KI.Raw,1)
    fig=figure;set(fig,'defaultAxesColorOrder',[[0.4 0 0.5]; [0 0 0]]);
    yyaxis left;
    plot(KI.Raw(iV,:),'b--d','MarkerEdgeColor','b','LineWidth',2); hold on;
    plot(KII.Raw(iV,:),'r--s','MarkerEdgeColor','r','LineWidth',2);
    plot(KIII.Raw(iV,:),'g--<','MarkerEdgeColor','g','LineWidth',2); hold off
    ylabel('K (MPa m^{0.5})');
    if min([KII.Avg(:); KI.Avg(:); KIII.Avg(:)])>0;     ylim([0 inf]);      end
    yyaxis right;
    plot(J.Raw(iV,:),'k--o','MarkerEdgeColor','k','LineWidth',2,'MarkerFaceColor','k');
    ylabel('J [J/m^2]');
    xlabel('Contour Number');
    legend(['K_{I} = ' num2str(KI.Ttrue(iV)) ' ± ' num2str(KI.Tdiv(iV) ) ' MPa\surdm' ],...
        ['K_{II} = ' num2str(KII.Ttrue(iV)) ' ± ' num2str(KII.Tdiv(iV) ) ' MPa\surdm' ],...
        ['K_{III} = ' num2str(KIII.Ttrue(iV)) ' ± ' num2str(KIII.Tdiv(iV) ) ' MPa\surdm' ],...
        ['J_{integral} = ' num2str(J.Ttrue(iV)) ' ± ' ...
        num2str(J.Tdiv(iV) ) ' J/m^2'],'location','northoutside','box','off');
    set(gcf,'WindowStyle','normal');box off;
    set(gcf,'position',[60,10,850,990]); xlim([0 length(KI.Avg)+1])
    saveas(gcf, [saveto '\' Ond '_' num2str(iV) '_KI, KII and J.fig']);
    saveas(gcf, [saveto '\' Ond '_' num2str(iV) '_KI, KII and J.tif']);    close
end

%% VALUES PER THICKNESS

fig=figure;set(fig,'defaultAxesColorOrder',[[0.4 0 0.5]; [0 0 0]]);
yyaxis left;
errorbar(thickness,KI.Ttrue,KI.Tdiv,'b--d','MarkerEdgeColor','b',...
    'LineWidth',2,'markersize',10); hold on;
errorbar(thickness,KII.Ttrue,KII.Tdiv,'r--s','MarkerEdgeColor','r',...
    'LineWidth',2,'markersize',10);
errorbar(thickness,KIII.Ttrue,KIII.Tdiv,'g--<','MarkerEdgeColor','g',...
    'LineWidth',2,'markersize',10); hold off
ylabel('K (MPa m^{0.5})');
if min([KII.Avg(:); KI.Avg(:); KIII.Avg(:)])>0;     ylim([0 inf]);      end
yyaxis right;
errorbar(thickness,J.Ttrue,J.Tdiv,'k--o','MarkerEdgeColor','k','LineWidth',2,...
    'MarkerFaceColor','k','markersize',10);
ylabel('J [J/m^2]');        %ylim([0 inf]);
xlabel(['Thickness (' unit ')']);
legend(['K_{I} = ' num2str(KI.true) ' ± ' num2str(KI.div) ' MPa\surdm' ],...
    ['K_{II} = ' num2str(KII.true) ' ± ' num2str(KII.div) ' MPa\surdm' ],...
    ['K_{III} = ' num2str(KIII.true) ' ± ' num2str(KIII.div) ' MPa\surdm' ],...
    ['J_{integral} = ' num2str(J.true) ' ± ' ...
    num2str(J.div) ' J/m^2'],'location','northoutside','box','off');
set(gcf,'WindowStyle','normal');box off;
set(gcf,'position',[60,10,850,990]);
saveas(gcf, [saveto '\' Ond '_Thickness_KI, KII and J.fig']);
saveas(gcf, [saveto '\' Ond '_Thickness_KI, KII and J.tif']);    close

% J
subplot(1,2,1)
errorbar(thickness,J.Ttrue,J.Tdiv,'k--o','MarkerEdgeColor','k',...
    'LineWidth',2,'MarkerFaceColor','k','markersize',10); hold on;
errorbar(thickness,J.K.Ttrue,J.K.Tdiv,'r--s','MarkerEdgeColor',...
    'r','LineWidth',2,'markersize',10); hold off
xlabel(['Thickness (' unit ')']);
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
    'location','northoutside','box','off');box off;
set(gcf,'WindowStyle','normal');    set(gcf,'position',[60 90 1298 882]);
saveas(gcf, [saveto '\' Ond '_Thickness_KI, KII and J.fig']);
saveas(gcf, [saveto '\' Ond '_Thickness_KI, KII and J.tif']);    close

save([saveto '\' Ond],'J','KI','KII','KIII')
%%
if ~exist(fullfile(saveto,[ok '_E12N.tif']),'file')
    fileid = fullfile(saveto, 'Abaqus_script.py');
    fileID = fopen(fileid,'w');
    ImportModules(fileID,saveto)
    SaveReport(fileID,saveto,ok)          %% save report
    printTiffs(fileID,saveto,'N',ok,3);   %% save tiff images for U and S
    printTiffs(fileID,saveto,'D',ok,3);   %% save tiff images for deformaed U and S
    fclose(fileID);
    PWD =pwd; cd(saveto);
    system(['abaqus cae ','noGUI','=Abaqus_script.py']); % Windows system
    cd(PWD);
end
