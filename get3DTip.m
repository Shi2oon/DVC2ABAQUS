function [xo,yo,dataum] = get3DTip(inData,xo,yo)
inData(:,1) = inData(:,1) - min(inData(:,1));
inData(:,2) = inData(:,2) - min(inData(:,2));
inData(:,3) = inData(:,3) - min(inData(:,3));

[ ~,dataum ] = reshapeData( inData );
U = (dataum.Uz(:,:,1));%.^2+dataum.Ux(:,:,1).^2+dataum.Uz(:,:,1).^2).^0.5;
imagesc(dataum.X1(1,:,1),dataum.Y1(:,1,1),U)%.^2+...
%         dataum.Uy(:,:,1).^2+dataum.Uz(:,:,1).^2).^0.5);
c=colorbar; c.Label.String = ['U_y'];
set(gca,'Ydir','normal');	axis image;colormap jet
title('Answer in the command line');
set(gcf,'WindowStyle','normal')
set(gcf,'position',[30 50 1300 950]);
title('Uy :: Select the Crack, start from crack tip');
if isempty(xo) || isempty(yo)
%{
opts.Interpreter = 'tex';       % Include the desired Default answer
opts.Default     = 'Y';         % Use the TeX interpreter to format the question
quest            = 'Is the crack is on your left?';
answer           = questdlg(quest,'Boundary Condition','Y','N', opts);
if  answer == 'N'
    inData(:,1) = abs(inData(:,1) - max(inData(:,1)));
    [ ~,dataum ] = reshapeData( inData );
    U = (dataum.Ux(:,:,1).^2+dataum.Uy(:,:,1).^2+dataum.Uz(:,:,1).^2).^0.5;
    imagesc(dataum.X1(1,:,1),dataum.Y1(:,1,1),U);
    c=colorbar; c.Label.String = ['Uy'];
    %     caxis([-5e-3 5e-3]);
    set(gca,'Ydir','normal');	axis image;colormap jet
    title('Answer in the command line');
    set(gcf,'WindowStyle','normal')
    set(gcf,'position',[30 50 1300 950]);
    title('Um :: Select the Crack, start from crack tip');
end
%}
[xo,yo] = ginput(2);

else
    xo = xo - min(inData(:,1));
    yo = yo - min(inData(:,2));
end
yo = [yo(1); yo(1)];

%% get excat from data in
xLin       = dataum.X1(1,:,1);
% [~, index] = min(abs(xLin-xo(1)));      xo(1) = xLin(index);
[~, index] = min(abs(xLin-xo(2)));      xo(2) = xLin(index);
% 
% yLin       = dataum.Y1(:,1,1);
% [~, index] = min(abs(yLin-yo(1)));      yo(1) = yLin(index);
% [~, index] = min(abs(yLin-yo(2)));      yo(2) = yLin(index);
line(xo,yo,'Color','w','LineStyle','-.','linewidth',2)

% if answer == 'Y'
    inData(:,1) = abs(inData(:,1) - max(inData(:,1)));
% else
%     inData(:,1) = inData(:,1) - min(inData(:,1));
% end
[ ~,dataum ] = reshapeData( inData );
dataum.X = dataum.X1;
dataum.Y = dataum.Y1;
dataum.Z = dataum.Z1;

%%
Zstep = unique(round(diff(unique(dataum.Z(:))),2));
Zstep(isnan(Zstep))=[];
if size(dataum.Z,3) == 2
    Zstep = min(Zstep)/2;
else
    Zstep = min(Zstep);
end
dataum.EdgeStepSize = min([diff(unique(dataum.X(:)))' diff(unique(dataum.Y(:)))' Zstep]);
dataum.MapStepSize = min([diff(unique(dataum.X(:)))' diff(unique(dataum.Y(:)))']);
if round(dataum.EdgeStepSize,3)==0
    dataum.EdgeStepSize = min([diff(unique(dataum.X(:)))' diff(unique(dataum.Y(:)))']);
end

title({'Crack Position (x,y) from tip to the end is at ',[ num2str(xo(1)) ' , ' ...
    num2str(yo(1)) ' and ' num2str(xo(2)) ' , ' num2str(yo(2))], ...
    ['Picture Frame is at ' ...
    num2str(min(unique(dataum.X1(:)))) ' , ' num2str(min(unique(dataum.Y1(:)))) ' and ' ...
    num2str(max(unique(dataum.X1(:)))) ' , ' num2str(max(unique(dataum.Y1(:))))],...
    [' Step size of ' num2str(dataum.MapStepSize)]});
end