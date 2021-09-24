function printDVC(saveto)
%%
[saveto,Ond] = fileparts(saveto);
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