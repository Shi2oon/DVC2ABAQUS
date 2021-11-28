function [ alldata,dataum ] = reshapeData( raw_data )
%PROCESS_DATA Summary of this function goes here
%   Detailed explanation goes here
% try
%     if size(raw_data,2) == 4
%         dataum.X1  = reshape(raw_data(:,1),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))));
%         dataum.Y1  = reshape(raw_data(:,2),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))));
%         dataum.Ux  = reshape(raw_data(:,3),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))));
%         dataum.Uy  = reshape(raw_data(:,4),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))));
%         alldata = [dataum.X1(:) dataum.Y1(:) dataum.Ux(:) dataum.Uy(:)];
%     else
%         dataum.X1  = reshape(raw_data(:,1),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))),length(unique(raw_data(:,3))));
%         dataum.Y1  = reshape(raw_data(:,2),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))),length(unique(raw_data(:,3))));
%         dataum.Z1  = reshape(raw_data(:,3),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))),length(unique(raw_data(:,3))));
%         dataum.Ux  = reshape(raw_data(:,4),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))),length(unique(raw_data(:,3))));
%         dataum.Uy  = reshape(raw_data(:,5),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))),length(unique(raw_data(:,3))));
%         dataum.Uz  = reshape(raw_data(:,6),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))),length(unique(raw_data(:,3))));
%         alldata = [dataum.X1(:) dataum.Y1(:) dataum.Z1(:) dataum.Ux(:) dataum.Uy(:) dataum.Uz(:)];
%     end
% catch
    
    x  = raw_data(:,1);
    y  = raw_data(:,2);
    ux = raw_data(:,3);
    uy = raw_data(:,4);
    
    xVec = unique(x);
    yVec = unique(y);
    
    % nDataPoints = length(x);
    % if length(xVec) <	length(raw_data(:,1))*0.5
    
    %Define grid
    [xMap,yMap] = meshgrid(xVec,yVec);
    [nRows, nCols] = size(xMap);
    
    % nGridPoints = length(xMap(:));
    
    uxMap = NaN(nRows, nCols); %Initialise
    uyMap = NaN(nRows, nCols); %Initialise
    
    if size(raw_data,2) == 6
        z  = raw_data(:,3);
        zVec = unique(z);
        ux = raw_data(:,4);
        uy = raw_data(:,5);
        uz = raw_data(:,6);
        [xMap,yMap,zMap] = meshgrid(xVec,yVec,zVec);
        [nRows, nCols , nDep] = size(xMap);
        uxMap = NaN(nRows, nCols, nDep); %Initialise
        uyMap = NaN(nRows, nCols, nDep); %Initialise
        uzMap = NaN(nRows, nCols, nDep); %Initialise
    end
    
    for iRow = 1:nRows % loop rows
        for iCol = 1:nCols % loop cols
            
            if size(raw_data,2) == 6 %% 3D
                for iDep = 1:nDep
                    xt = xMap(iRow,iCol,iDep);
                    yt = yMap(iRow,iCol,iDep);
                    zt = zMap(iRow,iCol,iDep);
                    idx = find(x==xt & y==yt & z==zt);
                    if ~isempty(idx)
                        uxt = ux(idx(1));
                        uyt = uy(idx(1));
                        uzt = uz(idx(1));
                        uxMap(iRow,iCol,iDep) = uxt;
                        uyMap(iRow,iCol,iDep) = uyt;
                        uzMap(iRow,iCol,iDep) = uzt;
                    end
                end
                
            else %% 2D
                xt = xMap(iRow,iCol);
                yt = yMap(iRow,iCol);
                idx = find(and(x==xt,y==yt)); %find linear index of point corresponding to xt,yt;
                if ~isempty(idx)
                    uxt = ux(idx(1));
                    uyt = uy(idx(1));
                    uxMap(iRow,iCol) = uxt;
                    uyMap(iRow,iCol) = uyt;
                end
            end
        end
    end
    
    dataum.X1 = xMap;
    dataum.Y1 = yMap;
    % dataum.Uy = uyMap;
    % dataum.Ux = uxMap;
    % threshold = 0.95;
    % [ uxMap ] = dispFieldSmoothing( uxMap, threshold );
    dataum.Ux = uxMap;
    % [ uyMap ] = dispFieldSmoothing( uyMap, threshold );
    dataum.Uy = uyMap;

alldata = [dataum.X1(:) dataum.Y1(:) dataum.Ux(:) dataum.Uy(:)];
if size(raw_data,2) == 6
    dataum.Z1 = zMap;
    dataum.Uz = uzMap;
    alldata = [dataum.X1(:) dataum.Y1(:) dataum.Z1(:) dataum.Ux(:) dataum.Uy(:) dataum.Uz(:)];
end
% end
%{
else
	disp('the data is highly non-uniform, I will now re-arrange the data');
%     scatter3(alldata(:,1), alldata(:,2), alldata(:,3),[],alldata(:,3)); view([0 90])
	Fx = scatteredInterpolant(raw_data(:,1), raw_data(:,2), raw_data(:,3),'natural','nearest');
	Fy = scatteredInterpolant(raw_data(:,1), raw_data(:,2), raw_data(:,4),'natural','nearest');
	X = linspace(min(raw_data(:,1)),max(raw_data(:,1)),300);
	Y = min(raw_data(:,2)):abs(X(2)-X(1)):max(raw_data(:,2));
	[dataum.X1,dataum.Y1] = meshgrid(X,Y);
	Ux = Fx(dataum.X1(:),dataum.Y1(:));
	Uy = Fy(dataum.X1(:),dataum.Y1(:));
    dataum.Ux = reshape(Ux,length(Y),length(X));
    dataum.Uy = reshape(Uy,length(Y),length(X));
	% scatter3(x(:),y(:),Ux,[],Ux); view([0 90])
end
%}


