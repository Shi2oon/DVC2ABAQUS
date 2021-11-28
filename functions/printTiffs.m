function printTiffs(fileID,folder,Stat,Unique,DIM)

fprintf(fileID,'session.viewports["Viewport: 1"].viewportAnnotationOptions.setValues(title=OFF, \n');
    fprintf(fileID,'    state=OFF, annotations=OFF, compass=OFF); \n');
    fprintf(fileID,'session.viewports["Viewport: 1"].viewportAnnotationOptions.setValues( \n');
    fprintf(fileID,'    legendFont="-*-verdana-medium-r-normal-*-*-120-*-*-p-*-*-*"); \n');
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.commonOptions.setValues(  \n');
    fprintf(fileID,'    visibleEdges=NONE);  \n');
    % U1
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(  \n');
    fprintf(fileID,'    variableLabel="U", outputPosition=NODAL, refinement=(COMPONENT, "U1"), ); \n');
    if Stat == 'N' % for normal view
        fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.display.setValues(plotState=(  \n');
        fprintf(fileID,'    CONTOURS_ON_UNDEF, ));  \n'); % not deform
    elseif Stat == 'D' % for deformed view
        fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.display.setValues(plotState=( \n');
        fprintf(fileID,'    CONTOURS_ON_DEF, ));\n'); % deformed view
    end
    SaveF = pythonFileName(fullfile(folder, [Unique '_U1' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');    
    % U2
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(  \n');
    fprintf(fileID,'    variableLabel="U", outputPosition=NODAL, refinement=(COMPONENT, "U2"), ); \n');
    SaveF = pythonFileName(fullfile(folder, [Unique '_U2' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
    %UM
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(  \n');
    fprintf(fileID,'    variableLabel="U", outputPosition=NODAL, refinement=(INVARIANT,"Magnitude"), ); \n');
    SaveF = pythonFileName(fullfile(folder, [Unique '_UM' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
    
    %U3
  if DIM==3
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(  \n');
    fprintf(fileID,'    variableLabel="U", outputPosition=NODAL, refinement=(COMPONENT, "U3"), ); \n');
    SaveF = pythonFileName(fullfile(folder, [Unique '_U3' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
  end
    
%% save stress   
    % Mises
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(\n');
    fprintf(fileID,'    variableLabel="S", outputPosition=INTEGRATION_POINT, refinement=(\n');
    fprintf(fileID,'    INVARIANT, "Mises"), ); \n');
    SaveF = pythonFileName(fullfile(folder, [Unique '_SM' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
    % S11
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(\n');
    fprintf(fileID,'    variableLabel="S", outputPosition=INTEGRATION_POINT, refinement=(\n');
    fprintf(fileID,'    COMPONENT, "S11"), );\n');
    SaveF = pythonFileName(fullfile(folder, [Unique '_S11' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
    % S12
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(\n');
    fprintf(fileID,'    variableLabel="S", outputPosition=INTEGRATION_POINT, refinement=(\n');
    fprintf(fileID,'    COMPONENT, "S12"), );\n');
    SaveF = pythonFileName(fullfile(folder, [Unique '_S12' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
    % S22
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(\n');
    fprintf(fileID,'    variableLabel="S", outputPosition=INTEGRATION_POINT, refinement=(\n');
    fprintf(fileID,'    COMPONENT, "S22"), );\n');
    SaveF = pythonFileName(fullfile(folder, [Unique '_S22' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
    
  if DIM==3
    % S13
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(\n');
    fprintf(fileID,'    variableLabel="S", outputPosition=INTEGRATION_POINT, refinement=(\n');
    fprintf(fileID,'    COMPONENT, "S13"), );\n');
    SaveF = pythonFileName(fullfile(folder, [Unique '_S13' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
    % S23
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(\n');
    fprintf(fileID,'    variableLabel="S", outputPosition=INTEGRATION_POINT, refinement=(\n');
    fprintf(fileID,'    COMPONENT, "S23"), );\n');
    SaveF = pythonFileName(fullfile(folder, [Unique '_S23' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
    % S33
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(\n');
    fprintf(fileID,'    variableLabel="S", outputPosition=INTEGRATION_POINT, refinement=(\n');
    fprintf(fileID,'    COMPONENT, "S33"), );\n');
    SaveF = pythonFileName(fullfile(folder, [Unique '_S33' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
  end
  
%% save stress   
    % E11
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(\n');
    fprintf(fileID,'    variableLabel="E", outputPosition=INTEGRATION_POINT, refinement=(\n');
    fprintf(fileID,'    COMPONENT, "E11"), );\n');
    SaveF = pythonFileName(fullfile(folder, [Unique '_E11' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
    % E12
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(\n');
    fprintf(fileID,'    variableLabel="E", outputPosition=INTEGRATION_POINT, refinement=(\n');
    fprintf(fileID,'    COMPONENT, "E12"), );\n');
    SaveF = pythonFileName(fullfile(folder, [Unique '_E12' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
    % E22
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(\n');
    fprintf(fileID,'    variableLabel="E", outputPosition=INTEGRATION_POINT, refinement=(\n');
    fprintf(fileID,'    COMPONENT, "E22"), );\n');
    SaveF = pythonFileName(fullfile(folder, [Unique '_E22' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
    
  if DIM==3
    % S13
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(\n');
    fprintf(fileID,'    variableLabel="E", outputPosition=INTEGRATION_POINT, refinement=(\n');
    fprintf(fileID,'    COMPONENT, "E13"), );\n');
    SaveF = pythonFileName(fullfile(folder, [Unique '_E13' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
    % S23
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(\n');
    fprintf(fileID,'    variableLabel="E", outputPosition=INTEGRATION_POINT, refinement=(\n');
    fprintf(fileID,'    COMPONENT, "E23"), );\n');
    SaveF = pythonFileName(fullfile(folder, [Unique '_E23' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
    % S33
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(\n');
    fprintf(fileID,'    variableLabel="E", outputPosition=INTEGRATION_POINT, refinement=(\n');
    fprintf(fileID,'    COMPONENT, "E33"), );\n');
    SaveF = pythonFileName(fullfile(folder, [Unique '_E33' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
  end