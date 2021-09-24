function [KI, KII,KIII, JK] = read3DDATAbaqus(filrname)
    
fid = fopen(filrname,'rt') ;
S = textscan(fid,'%s','Delimiter','\n');
S = S{1} ;

%% J 
idxS = strfind(S, 'INCREMENT    20 SUMMARY');
idx1 = find(not(cellfun('isempty', idxS)));
idx1 = idx1(end);
idxS = strfind(S, 'LABELS REFERENCED IN THE ABOVE TABLE');
idx2 = find(not(cellfun('isempty', idxS)));
idx2 = idx2(end);

%% KII and KII KIII
K_in = S(idx1+18:idx2-5);
for iv = length(K_in):-1:1
    if isempty(K_in{iv,1})
        K_in(iv)=[];
    end
end

B  = convertCharsToStrings(K_in);
count=0;
for iv = 1:length(B)
    C{iv}  = textscan(B(iv),'%s');
    C{iv}  = C{iv}{1,1};
    format = 0;
    count = count+1;
    for ic = 1:length(C{iv})
        Conv = str2double(C{iv}(ic));
        if isnan(Conv)==0
            format = format+1;
            OutKJ(count,format) = str2double(C{iv}(ic));
        end
    end
end
count=1;
for iv = 2:4:size(OutKJ,1)
    KI(count:count+size(OutKJ,2)-1)  = OutKJ(iv,:);          
    KII(count:count+size(OutKJ,2)-1) = OutKJ(iv+1,:);    
    KIII(count:count+size(OutKJ,2)-1) = OutKJ(iv+2,:);  
    JK(count:count+size(OutKJ,2)-1)  = OutKJ(iv+3,:);    
    count = count+size(OutKJ,2);
end

%% correct
KI(KI==0)  =[];            KII(KII==0)=[];              KIII(KIII==0)=[];                     
JK(JK==0)  =[];            JK  = abs(JK(:));
% if mean(sign(KI)) ~= mean(sign(KII))
%     if mean(abs(KI)) > mean(abs(KII))
%         KI = abs(KI);       
%     else
%         KII = abs(KII);      
%     end
% else
%     KI  = abs(KI(:));       KII = abs(KII(:));
% end
fclose('all');

close all;
set(0,'defaultAxesFontSize',25);       set(0,'DefaultLineMarkerSize',14) 
plot(abs(KI),'r.','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r');     hold on
plot(abs(KII),'b.','MarkerEdgeColor','b','LineWidth',1.5,'MarkerFaceColor','b'); 
plot(abs(KIII),'k.','MarkerEdgeColor','k','LineWidth',1.5,'MarkerFaceColor','k');    hold off;
xlabel('Slices');       ylabel('K (MPa\surdm)')
title ([ 'K_I = ' num2str(round(trimmean(abs(KI),50),3)) ' MPa\surdm, K_{II} = '...
    round(num2str(trimmean(abs(KII),50),3))...
    ' MPa\surdm, K_{III} = '  round(num2str(trimmean(abs(KIII),50),3)) ' MPa\surdm']);
legend('K_{I}','K_{II}','K_{III}','location','best'); 
set(gcf,'position',[600,100,950,650]);   xl = xlim;  xlim([0 xl(2)+1])
saveas(gcf, [erase(filrname,'.dat') '_K.fig']); 
saveas(gcf, [erase(filrname,'.dat') '_K.png']); close

plot(abs(JK),'r.','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r');    
xlabel('Slices');       ylabel('J (J/m^2)')
set(gcf,'position',[600,100,950,650]);   xl = xlim;  xlim([0 xl(2)+1])
saveas(gcf, [erase(filrname,'.dat') '_J.fig']); 
saveas(gcf, [erase(filrname,'.dat') '_J.png']); close