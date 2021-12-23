function  writeDataToFile(requestedComponents,requestWriteAll,printHeader,inDir,outFolder,data_col,dataSize,cycleName,operations,eulerAngles)
% Write columnwise 3D (2D) data to file, as [coordinates data] as if a .DAT
% file. A 3-line header containing processing details is optional.
%
% (c) MSL Jordan, University of Oxford 2015
%% 1. Definitions
% disp('*********** WARNING: folder Error ***********')
outDir = [inDir outFolder '\'];
c=clock;
component = {'U_x','U_y','U_z','|U|'};

switch size(dataSize,2)
    case 3
        dims=3;
        inCo = {1:3};   %Index for coordinate columns
        
        %Header strings
        string2 = ['|       x/mm|        y/mm|        z/mm|      ' component{1} '/mm|      ' component{2} '/mm|      ' component{3} '/mm|'];
        string3 = ['Properties: Size /voxel I=' sprintf('%03d',dataSize(1)) ', J=' sprintf('%03d',dataSize(2)) ', K=' sprintf('%03d',dataSize(3))]; 
    case 2
        dims=2;
        inCo = {1:2};   %Index for coordinate columns
        
        %Header strings
        string2 = ['|       x/mm|        y/mm|      ' component{1} '/mm|      ' component{2} '/mm|'];
        string3 = ['Properties: Size /voxel I=' sprintf('%03d',dataSize(1)) ', J=' sprintf('%03d',dataSize(2))]; 
end

%% 2. Create output folder
if ~exist(outDir,'dir')
    mkdir(inDir,outFolder);
    disp(['Saving rotated data in ' outDir])
else
    disp('Output folder already exists; data will be overwritten.')
    disp('=>Hit any key to continue')
%     pause
    disp(['Saving rotated data in ' outDir])
end

%% 3. Write Data
if requestWriteAll
    disp('Writing DVC-all');
    filename = [cycleName '-All.txt'];
    outPath = [outDir filename];
    
    string1 = ['Dataset: ' filename ' rotated data calculated at ' num2str(c(4)) ':' sprintf('%02d',c(5)) ' on ' num2str(c(3)) '/' num2str(c(2)) '/' num2str(c(1)) ', Operations: [' num2str(operations) '], Euler_angles = ' num2str(radtodeg(eulerAngles),'% 1.3f ')];
    
    fileID3 = fopen(outPath,'w');
    if printHeader
        fprintf(fileID3,'%s\r\n',string1);
        fprintf(fileID3,'%s\r\n',string2);
        fprintf(fileID3,'%s\r\n',string3);
    end
    
    switch dims
        case 3
            fprintf(fileID3,'%11.6E\t%11.6E\t%11.6E\t%11.6E\t%11.6E\t%11.6E \r\n',data_col');
        case 2
            fprintf(fileID3,'%11.6E\t%11.6E\t%11.6E\t%11.6E \r\n',data_col');
    end
    fclose(fileID3);
end

if numel(requestedComponents)>0
    disp('Writing displacement components to file');
    
    % For compatibility with "Displacement_difference.exe" (Delphi file)
    % require 4 column format:
    %
    % string1//
    % string2//
    % string3//
    % x|y|z|Ui

    writeData = zeros(size(data_col,1),4);
    writeData(:,inCo) = data_col(:,inCo);

    for i = requestedComponents
        filename = [cycleName '-' component{i} '.txt'];
        outPath = [outDir filename];

        string1 = [cycleName '-' component{i} ' rotated data calculated at ' num2str(c(4)) ':' sprintf('%02d',c(5)) ' on ' num2str(c(3)) '/' num2str(c(2)) '/' num2str(c(1)) ', Operations: ' operations ', Euler_angles = ' num2str(radtodeg(t'),'% 1.3f ')];

%         writeData(:,inCo) = data_col(:,inCo);
        writeData(:,4) = data_col(:,i+dims);

        fileID3 = fopen(outPath,'w');
        if printHeader
            fprintf(fileID3,'%s\r\n',string1);       %probably change to one line
            fprintf(fileID3,'%s\r\n',string2);
            fprintf(fileID3,'%s\r\n',string3);
        end
        fprintf(fileID3,'%11.6E\t%11.6E\t%11.6E\t%11.6E \r\n',writeData');
        fclose(fileID3);
    end
end

disp('Data written to file.')
end
