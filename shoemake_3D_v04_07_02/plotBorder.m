function plotBorder(data_col, mapID)
%Plots a bounding box around data on a 2D orthoslice
%
% (c) MSL Jordan, University of Oxford, 2015

%% Definitions

dims = size(data_col,2);

switch dims
    case {2 4}
        X_temp = 1;
        Y_temp = 2;
    case {3 6}
        switch mapID(1)
            case 1
                X_temp = 2;
                Y_temp = 3;
            case 2
                X_temp = 1;
                Y_temp = 3;        
            case 3
                X_temp = 1;
                Y_temp = 2;  
            otherwise
                disp '##ERROR: map direction (mapID(1) invalid'; ERROR
        end
    otherwise
        disp '##ERROR: data_col invalid number of columns'; ERROR
end

Hi = max(data_col(:,1:dims));
Lo = min(data_col(:,1:dims));  

%% Compose border coordinates
border = [...
    Hi(X_temp) Hi(Y_temp);...
    Hi(X_temp) Lo(Y_temp);...
    Lo(X_temp) Lo(Y_temp);...
    Lo(X_temp) Hi(Y_temp);...
    Hi(X_temp) Hi(Y_temp)];  

%% Add to current plot
hold on
plot(border(:,1),border(:,2),'c-');
