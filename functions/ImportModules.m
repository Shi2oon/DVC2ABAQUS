function ImportModules(fileID,folder)
    fprintf(fileID,'from __future__ import division  \n');
    % ABAQUS
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
    fprintf(fileID,'\n');
    % NOT ABAQUS BUT EXISTING LIBS 
    fprintf(fileID,'from operator import itemgetter \n');
    fprintf(fileID,'from math import ceil \n');
    fprintf(fileID,'import numpy  as np \n');
    fprintf(fileID,'import os \n');
    fprintf(fileID,'\n');
    
   % set directory to output file
   FoldrOut = pythonFileName(folder);
   fprintf(fileID,'os.chdir(r"%s");\n',FoldrOut);
end