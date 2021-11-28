function SaveReport(fileID,folder,MatPunique)
    %% save report
    Saveodb = pythonFileName(fullfile(folder, [MatPunique '.odb']));
    fprintf(fileID,'odb = session.openOdb("%s"); \n',Saveodb);
    fprintf(fileID,'session.viewports["Viewport: 1"].setValues(displayedObject=odb);\n');
    fprintf(fileID,'session.fieldReportOptions.setValues(sort=ASCENDING);\n');
    fprintf(fileID,'session.fieldReportOptions.setValues(columnLayout=SEPARATE_TABLES);\n');
    fprintf(fileID,'session.writeFieldReport( \n');
    SaveReport = pythonFileName(fullfile(folder, [MatPunique '.rpt']));
    fprintf(fileID,'    fileName="%s", \n',SaveReport);
	fprintf(fileID,'    append=ON, sortItem="Node Label", odb=odb, step=0, frame=1, \n');
	fprintf(fileID,'    outputPosition=ELEMENT_NODAL, variable=(("S", INTEGRATION_POINT, ((\n');
	fprintf(fileID,'    INVARIANT, "Mises"), (COMPONENT, "S11"), (COMPONENT, "S22"), (  \n');
	fprintf(fileID,'    COMPONENT, "S12"), )), ));\n');

end