function importDat(fileToRead)
%IMPORTDAT(FILETOREAD)
%  Imports data from the specified file
%  FILETOREAD:  file to read
%
% Based on MATLAB auto-generated script. Modified by MSL Jordan
% (c) MSL Jordan, University of Oxford 2015


% Import the file
newData = importdata(fileToRead);

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData);
% size(vars,1)

for i = [1 size(vars,1)]
    assignin('caller', vars{i}, newData.(vars{i})); %%N.B. shoemake works using 'base'
end

if size(vars,1)==2
    assignin('caller', 'colheaders', strsplit(newData.textdata{end,:})); %%N.B. shoemake works using 'base'
end 
