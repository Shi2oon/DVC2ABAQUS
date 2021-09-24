function Greate3DCracked_v2(folder,Len,A,offset)
% Len = [min(X) max(X); min(y) max(Y); z cracklength];
PWD = pwd;
cd(folder);
fileid = [folder, '\Abaqus_script_1.py'];
fileID = fopen(fileid,'w');
fprintf(fileID,'from __future__ import division  \n');
fprintf(fileID,'from abaqus import * \n');
fprintf(fileID,'from abaqusConstants import * \n');
fprintf(fileID,'from part import * \n');
fprintf(fileID,'from material import * \n');
fprintf(fileID,'from section import * \n');
fprintf(fileID,'from assembly import * \n');
fprintf(fileID,'from step import * \n');
fprintf(fileID,'from interaction import * \n');
fprintf(fileID,'from load import * \n');
fprintf(fileID,'from mesh import * \n');
fprintf(fileID,'from optimization import * \n');
fprintf(fileID,'from job import * \n');
fprintf(fileID,'from sketch import * \n');
fprintf(fileID,'from visualization import * \n');
fprintf(fileID,'from connectorBehavior import * \n');
fprintf(fileID,'from odbAccess import * \n');
fprintf(fileID,'import regionToolset\n');
fprintf(fileID,'import displayGroupMdbToolset as dgm\n');
fprintf(fileID,'import xyPlot\n');
fprintf(fileID,'import displayGroupOdbToolset as dgo\n');
fprintf(fileID,'import connectorBehavior\n');
fprintf(fileID,'from operator import itemgetter \n');
fprintf(fileID,'from math import ceil \n');
fprintf(fileID,'import numpy  as np \n');
fprintf(fileID,'import os \n');
FoldrOut = pythonFileName(folder);
fprintf(fileID,'os.chdir(r"%s");\n',FoldrOut);
% sketch
fprintf(fileID,'mdb.models["Model-1"].ConstrainedSketch(name="__profile__", sheetSize=200.0)\n');
if Len(4,1) > Len(4,2) % crack on ur left
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].Line(point1=(%f, %f), point2=(\n',Len(1,1),Len(2,2));
    fprintf(fileID,'    %f, %f))\n',Len(1,2),Len(2,2));
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].HorizontalConstraint(\n');
    fprintf(fileID,'    addUndoState=False, entity=\n');
    fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[2])\n');
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].Line(point1=(%f, %f), point2=\n',Len(1,2),Len(2,2));
    fprintf(fileID,'    (%f, %f))\n',Len(1,2),Len(2,1));
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].VerticalConstraint(addUndoState=\n');
    fprintf(fileID,'    False, entity=mdb.models["Model-1"].sketches["__profile__"].geometry[3])\n');
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].PerpendicularConstraint(\n');
    fprintf(fileID,'    addUndoState=False, entity1=\n');
    fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[2], entity2=\n');
    fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[3])\n');
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].Line(point1=(%f, %f), point2=(\n',Len(1,2),Len(2,1));
    fprintf(fileID,'    %f, %f))\n',Len(1,1),Len(2,1));
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].HorizontalConstraint(\n');
    fprintf(fileID,'    addUndoState=False, entity=\n');
    fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[4])\n');
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].PerpendicularConstraint(\n');
    fprintf(fileID,'    addUndoState=False, entity1=\n');
    fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[3], entity2=\n');
    fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[4])\n');
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].Line(point1=(%f, %f), point2=(\n',Len(1,1),Len(2,1));
    fprintf(fileID,'    %f, %f))\n',Len(1,1),Len(5,2)-Len(6,1)/2);
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].VerticalConstraint(addUndoState=\n');
    fprintf(fileID,'    False, entity=mdb.models["Model-1"].sketches["__profile__"].geometry[5])\n');
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].PerpendicularConstraint(\n');
    fprintf(fileID,'    addUndoState=False, entity1=\n');
    fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[4], entity2=\n');
    fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[5])\n');
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].Line(point1=(%f, %f), point2=(\n',Len(1,1),Len(5,2)-Len(6,1)/2);
    fprintf(fileID,'    %f, %f))\n',Len(4,1),Len(5,2));
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].Line(point1=(%f, %f), point2=\n',Len(4,1),Len(5,2));
    fprintf(fileID,'    (%f, %f))\n',Len(1,1),Len(5,2)+Len(6,1)/2);
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].Line(point1=(%f, %f), point2=(\n',Len(1,1),Len(5,2)+Len(6,1)/2);
    fprintf(fileID,'    %f, %f))\n',Len(1,1),Len(2,2));
else % crack on ur right
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].Line(point1=(%f, %f), point2=(\n',Len(1,2),Len(2,2));
    fprintf(fileID,'    %f, %f))\n',Len(1,1),Len(2,2));
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].HorizontalConstraint(\n');
    fprintf(fileID,'    addUndoState=False, entity=\n');
    fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[2])\n');
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].Line(point1=(%f, %f), point2=\n',Len(1,1),Len(2,2));
    fprintf(fileID,'    (%f, %f))\n',Len(1,1),Len(2,1));
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].VerticalConstraint(addUndoState=\n');
    fprintf(fileID,'    False, entity=mdb.models["Model-1"].sketches["__profile__"].geometry[3])\n');
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].PerpendicularConstraint(\n');
    fprintf(fileID,'    addUndoState=False, entity1=\n');
    fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[2], entity2=\n');
    fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[3])\n');
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].Line(point1=(%f, %f), point2=(\n',Len(1,1),Len(2,1));
    fprintf(fileID,'    %f, %f))\n',Len(1,2),Len(2,1));
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].HorizontalConstraint(\n');
    fprintf(fileID,'    addUndoState=False, entity=\n');
    fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[4])\n');
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].PerpendicularConstraint(\n');
    fprintf(fileID,'    addUndoState=False, entity1=\n');
    fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[3], entity2=\n');
    fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[4])\n');
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].Line(point1=(%f, %f), point2=(\n',Len(1,2),Len(2,1));
    fprintf(fileID,'    %f, %f))\n',Len(1,2),Len(5,2)-Len(6,1)/2);
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].VerticalConstraint(addUndoState=\n');
    fprintf(fileID,'    False, entity=mdb.models["Model-1"].sketches["__profile__"].geometry[5])\n');
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].PerpendicularConstraint(\n');
    fprintf(fileID,'    addUndoState=False, entity1=\n');
    fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[4], entity2=\n');
    fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[5])\n');
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].Line(point1=(%f, %f), point2=(\n',Len(1,2),Len(5,2)-Len(6,1)/2);
    fprintf(fileID,'    %f, %f))\n',Len(4,1),Len(5,2));
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].Line(point1=(%f, %f), point2=\n',Len(4,1),Len(5,2));
    fprintf(fileID,'    (%f, %f))\n',Len(1,2),Len(5,2)+Len(6,1)/2);
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].Line(point1=(%f, %f), point2=(\n',Len(1,2),Len(5,2)+Len(6,1)/2);
    fprintf(fileID,'    %f, %f))\n',Len(1,2),Len(2,2));
end

fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].VerticalConstraint(addUndoState=\n');
fprintf(fileID,'    False, entity=mdb.models["Model-1"].sketches["__profile__"].geometry[8])\n');
fprintf(fileID,'mdb.models["Model-1"].Part(dimensionality=THREE_D, name="Part-1", type=\n');
fprintf(fileID,'    DEFORMABLE_BODY)\n');
fprintf(fileID,'mdb.models["Model-1"].parts["Part-1"].BaseSolidExtrude(depth=%f, sketch=\n',Len(3,1));
fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"])\n');
fprintf(fileID,'del mdb.models["Model-1"].sketches["__profile__"]\n');

%{
fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].HorizontalConstraint(\n');
fprintf(fileID,'    addUndoState=False, entity=\n');
fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[6])\n');
fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].PerpendicularConstraint(\n');
fprintf(fileID,'    addUndoState=False, entity1=\n');
fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[2], entity2=\n');
fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[6])\n');
fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].CoincidentConstraint(\n');
fprintf(fileID,'    addUndoState=False, entity1=\n');
fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].vertices[4], entity2=\n');
fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].geometry[2])\n');
fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].EqualDistanceConstraint(\n');
fprintf(fileID,'    addUndoState=False, entity1=\n');
fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].vertices[0], entity2=\n');
fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].vertices[1], midpoint=\n');
fprintf(fileID,'    mdb.models["Model-1"].sketches["__profile__"].vertices[4])\n');
fprintf(fileID,'mdb.models["Model-1"].parts["Part-1"].PartitionFaceBySketch(faces=\n');
fprintf(fileID,'    mdb.models["Model-1"].parts["Part-1"].faces.getSequenceFromMask(("[#10 ]", \n');
fprintf(fileID,'    ), ), sketch=mdb.models["Model-1"].sketches["__profile__"], sketchUpEdge=\n');
fprintf(fileID,'    mdb.models["Model-1"].parts["Part-1"].edges[0])\n');
fprintf(fileID,'del mdb.models["Model-1"].sketches["__profile__"]\n');
fprintf(fileID,'mdb.models["Model-1"].parts["Part-1"].PartitionCellBySweepEdge(cells=\n');
fprintf(fileID,'    mdb.models["Model-1"].parts["Part-1"].cells.getSequenceFromMask(("[#1 ]", \n');
fprintf(fileID,'    ), ), edges=(mdb.models["Model-1"].parts["Part-1"].edges[13], ), sweepPath=\n');
fprintf(fileID,'    mdb.models["Model-1"].parts["Part-1"].edges[2])\n');
%}

fprintf(fileID,'mdb.models["Model-1"].Material(name="Cracked")\n');
if A.type == 'A'
    A.C = A.Stiffness*offset^2;
    fprintf(fileID,'mdb.models["Model-1"].materials["Cracked"].Elastic(table=((%f, %f, \n',A.C(1,1),A.C(1,2));
    fprintf(fileID,'    %f, %f, %f, %f, %f, %f, %f, %f, %f, \n',...
        A.C(1,3),A.C(1,4),A.C(1,5),A.C(1,6),A.C(2,2),A.C(2,3),A.C(2,4),A.C(2,5),A.C(2,6));
    fprintf(fileID,'    %f, %f, %f, %f, %f, %f, %f, %f, %f, \n',...
        A.C(3,3),A.C(3,4),A.C(3,5),A.C(3,6),A.C(4,4),A.C(4,5),A.C(4,6),A.C(5,5),A.C(5,6));
    fprintf(fileID,'    %f), ), type=ANISOTROPIC)\n',A.C(6,6));
elseif A.type == 'E'
    A.E = A.E*offset^2;
    fprintf(fileID,'mdb.models["Model-1"].materials["Cracked"].Elastic(table=((%f, %f),))\n',A.E,A.nu);
elseif A.type == 'R'
    A.E = A.E*offset^2;       A.yield = A.yield*offset^2;
    fprintf(fileID,'	mdb.models["Model-1"].materials["Cracked"].DeformationPlasticity(table=((%f,  \n', A.E);
    fprintf(fileID,'		%f, %f, %f, %f), )) \n',A.nu,A.yield,A.Exponent, A.Yield_offset);
end

fprintf(fileID,'mdb.models["Model-1"].HomogeneousSolidSection(material="Cracked", name=\n');
fprintf(fileID,'    "Section-1", thickness=None)\n');
fprintf(fileID,'mdb.models["Model-1"].parts["Part-1"].Set(cells=\n');
fprintf(fileID,'    mdb.models["Model-1"].parts["Part-1"].cells.getSequenceFromMask(("[#1 ]", \n');
fprintf(fileID,'    ), ), name="Set-1")\n');
fprintf(fileID,'mdb.models["Model-1"].parts["Part-1"].SectionAssignment(offset=0.0, \n');
fprintf(fileID,'    offsetField="", offsetType=MIDDLE_SURFACE, region=\n');
fprintf(fileID,'    mdb.models["Model-1"].parts["Part-1"].sets["Set-1"], sectionName=\n');
fprintf(fileID,'    "Section-1", thicknessAssignment=FROM_SECTION)\n');
fprintf(fileID,'mdb.models["Model-1"].rootAssembly.DatumCsysByDefault(CARTESIAN)\n');
fprintf(fileID,'mdb.models["Model-1"].rootAssembly.Instance(dependent=OFF, name="Part-1-1", \n');
fprintf(fileID,'    part=mdb.models["Model-1"].parts["Part-1"])\n');
fprintf(fileID,'mdb.models["Model-1"].parts["Part-1"].MaterialOrientation(\n');
fprintf(fileID,'    additionalRotationType=ROTATION_NONE, axis=AXIS_1, fieldName="", localCsys=\n');
fprintf(fileID,'    None, orientationType=GLOBAL, region=Region(\n');
fprintf(fileID,'    cells=mdb.models["Model-1"].parts["Part-1"].cells.getSequenceFromMask(\n');
fprintf(fileID,'    mask=("[#1 ]", ), )), stackDirection=STACK_3)\n');
fprintf(fileID,'mdb.models["Model-1"].rootAssembly.regenerate()\n');
fprintf(fileID,'mdb.models["Model-1"].StaticStep(name="Step-1", previous="Initial")\n');
fprintf(fileID,'mdb.models["Model-1"].rootAssembly.Set(faces=\n');
fprintf(fileID,'    mdb.models["Model-1"].rootAssembly.instances["Part-1-1"].faces.getSequenceFromMask(\n');
fprintf(fileID,'    ("[#1 ]", ), ), name="Set-1")\n');
fprintf(fileID,'mdb.models["Model-1"].rootAssembly.engineeringFeatures.assignSeam(regions=\n');
fprintf(fileID,'    mdb.models["Model-1"].rootAssembly.sets["Set-1"])\n');

%% Meshing
% fprintf(fileID,'mdb.models["Model-1"].rootAssembly.seedPartInstance(deviationFactor=0.1, \n');
% fprintf(fileID,'    minSizeFactor=0.1, regions=(\n');
% fprintf(fileID,'    mdb.models["Model-1"].rootAssembly.instances["Part-1-1"], ), size=%f)\n',Len(6,1));
fprintf(fileID,'mdb.models["Model-1"].rootAssembly.seedPartInstance(deviationFactor=0.1, \n');
fprintf(fileID,'    minSizeFactor=0.1, regions=(\n');
fprintf(fileID,'    mdb.models["Model-1"].rootAssembly.instances["Part-1-1"], ), size=%f)\n',Len(3,2));
fprintf(fileID,'mdb.models["Model-1"].rootAssembly.seedEdgeBySize(constraint=FINER, \n');
fprintf(fileID,'    deviationFactor=0.1, edges=\n');
fprintf(fileID,'    mdb.models["Model-1"].rootAssembly.instances["Part-1-1"].edges.getSequenceFromMask(\n');
fprintf(fileID,'    ("[#20000 ]", ), ), minSizeFactor=0.1, size=%f)\n',Len(6,1));
% generate mesh
fprintf(fileID,'mdb.models["Model-1"].rootAssembly.generateMesh(regions=(\n');
fprintf(fileID,'    mdb.models["Model-1"].rootAssembly.instances["Part-1-1"], ))\n');

% Job
fprintf(fileID,'mdb.Job(atTime=None, contactPrint=OFF, description="", echoPrint=OFF, \n');
fprintf(fileID,'    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, \n');
fprintf(fileID,'    memory=90, memoryUnits=PERCENTAGE, model="Model-1", modelPrint=OFF, \n');
fprintf(fileID,'    multiprocessingMode=DEFAULT, name="CrackedModel", nodalOutputPrecision=SINGLE, \n');
fprintf(fileID,'    numCpus=8, numDomains=8, numGPUs=8, queue=None, resultsFormat=ODB, scratch=\n');
fprintf(fileID,'    "", type=ANALYSIS, userSubroutine="", waitHours=0, waitMinutes=0)\n');
% fprintf(fileID,'mdb.jobs["CrackedModel"].submit(consistencyChecking=OFF)\n');

% define the crack
if Len(4,1) > Len(4,2) % crack on ur left
    fprintf(fileID,'mdb.models["Model-1"].rootAssembly.engineeringFeatures.ContourIntegral(\n');
    fprintf(fileID,'    collapsedElementAtTip=NONE, crackFront=Region(\n');
    fprintf(fileID,'    edges=mdb.models["Model-1"].rootAssembly.instances["Part-1-1"].edges.getSequenceFromMask(\n');
    fprintf(fileID,'    mask=("[#4000 ]", ), )), crackTip=Region(\n');
    fprintf(fileID,'    edges=mdb.models["Model-1"].rootAssembly.instances["Part-1-1"].edges.getSequenceFromMask(\n');
    fprintf(fileID,'    mask=("[#4000 ]", ), )), extensionDirectionMethod=Q_VECTORS, \n');
    fprintf(fileID,'    midNodePosition=0.5, name="Crack-1", qVectors=((\n');
    fprintf(fileID,'    mdb.models["Model-1"].rootAssembly.instances["Part-1-1"].vertices[6], \n');
    fprintf(fileID,'    mdb.models["Model-1"].rootAssembly.instances["Part-1-1"].InterestingPoint(\n');
    fprintf(fileID,'    mdb.models["Model-1"].rootAssembly.instances["Part-1-1"].edges[7], \n');
    fprintf(fileID,'    MIDDLE)), ), symmetric=OFF)\n');
else
    fprintf(fileID,'mdb.models["Model-1"].rootAssembly.engineeringFeatures.ContourIntegral(\n');
    fprintf(fileID,'    collapsedElementAtTip=NONE, crackFront=Region(\n');
    fprintf(fileID,'    edges=mdb.models["Model-1"].rootAssembly.instances["Part-1-1"].edges.getSequenceFromMask(\n');
    fprintf(fileID,'    mask=("[#2 ]", ), )), crackTip=Region(\n');
    fprintf(fileID,'    edges=mdb.models["Model-1"].rootAssembly.instances["Part-1-1"].edges.getSequenceFromMask(\n');
    fprintf(fileID,'    mask=("[#2 ]", ), )), extensionDirectionMethod=Q_VECTORS, midNodePosition=\n');
    fprintf(fileID,'    0.5, name="Crack-1", qVectors=((\n');
    fprintf(fileID,'    mdb.models["Model-1"].rootAssembly.instances["Part-1-1"].vertices[12], \n');
    fprintf(fileID,'    mdb.models["Model-1"].rootAssembly.instances["Part-1-1"].InterestingPoint(\n');
    fprintf(fileID,'    mdb.models["Model-1"].rootAssembly.instances["Part-1-1"].edges[16], \n');
    fprintf(fileID,'    MIDDLE)), ), symmetric=OFF)\n');
end

fprintf(fileID,'mdb.models["Model-1"].HistoryOutputRequest(contourIntegral="Crack-1",\n');
fprintf(fileID,'    createStepName="Step-1", name="H-Output-2", numberOfContours=%d, rebar=\n',Len(6,2));
fprintf(fileID,'    EXCLUDE, sectionPoints=DEFAULT)\n');
fprintf(fileID,'mdb.models["Model-1"].HistoryOutputRequest(contourIntegral="Crack-1", \n');
fprintf(fileID,'    contourType=K_FACTORS, createStepName="Step-1", name="H-Output-3", \n');
fprintf(fileID,'    numberOfContours=%d, rebar=EXCLUDE, sectionPoints=DEFAULT)\n',Len(6,2));
fprintf(fileID,'mdb.jobs["CrackedModel"].submit(consistencyChecking=OFF)\n');
fclose(fileID);

system(['abaqus cae ','noGUI','=Abaqus_script_1.py']); % Windows system
cd(PWD)
delete([folder '\CrackedModel.log']);
delete([folder '\CrackedModel.msg']);
delete([folder '\CrackedModel.odb']);
delete([folder '\CrackedModel.prt']);
delete([folder '\CrackedModel.sim']);
delete([folder '\CrackedModel.sta']);
delete([folder '\CrackedModel.com']);