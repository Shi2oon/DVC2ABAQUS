function [RMS, Ave] = stdevAndAverage(cleanedData,dataSize,z_range,dx_mm)
% (c) MSL Jordan, University of Oxford 2015

%{ 
%TEST VARS
cleanedData = big_new_data(:,4:6);
% data_name = 'U_y';
% size = [lx ly lz];
z_range = [10 lz-10];

% Z_coord feature not implemented
%}

z_length = z_range(2)-z_range(1)+1;

% RMS and Ave have column headings: 
%
% ||Slice_number | Z_coord|Blank!|<Ux>|<Uy>|<Uz>|<U>|{SD(U)/<U>|}|

nCols = size(cleanedData,2);
% RMS = zeros(z_length,nCols);
% Ave=RMS;

RMS = zeros(z_length,3+nCols+2);
RMS(:,1) = (z_range(1):z_range(2));
RMS(:,2) = RMS(:,1)*dx_mm;
Ave=RMS;

for j = 1:(nCols)
    cData_temp=reshape(cleanedData(:,j),dataSize);
    
    for i = z_range(1):z_range(2)
        entry = i-z_range(1)+1;

        %slice_range = [((i-1)*slice_size)+1 i*slice_size];
    %     RMS(entry) = nanstd(cleaned_data(:,i));           %for 2-D work
    %     Ave(entry) = nanmean(cleaned_data(:,i));
        RMS(entry,j+3) = nanstd(nanstd(cData_temp(:,:,i)));   %for 3-D work
        Ave(entry,j+3) = nanmean(nanmean(cData_temp(:,:,i)));
    end
end

RMS(:,7) = sqrt(RMS(:,4).^2+RMS(:,5).^2+RMS(:,6).^2);
RMS(:,8) = sqrt((RMS(:,4)./Ave(:,4)).^2+(RMS(:,5)./Ave(:,5)).^2+(RMS(:,6)./Ave(:,6)).^2);
Ave(:,7) = sqrt(Ave(:,4).^2+Ave(:,5).^2+Ave(:,6).^2);
%% Comparing the average values (in pixels)
%{
disp 'averages'
single_average = nanmean(nanmean(nanmean(cleaned_data)))/1.8e-3
layer_average = nanmean(Ave)/1.8e-3
%}
end
