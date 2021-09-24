function [Nodes,Elements] = readNodes(filrname)
%% iinp file
fid = fopen(filrname,'rt') ;
AllText = textscan(fid,'%s','Delimiter','\n');
AllText = AllText{1} ;fclose('all');
%% Get Nodes and elements
% Nodes: number x y
idxS = strfind(AllText, 'Node');
idx1 = find(not(cellfun('isempty', idxS)));
idxS = strfind(AllText, 'Element');
idx2 = find(not(cellfun('isempty', idxS)));
idxS = strfind(AllText, 'End');
idx3 = find(not(cellfun('isempty', idxS)));
idxS = strfind(AllText, 'Nset');
idx4 = find(not(cellfun('isempty', idxS)));
if ~isempty(idx4); idx3(1) = idx4(1); end

Nodes = AllText(idx1+1:idx2-1) ;
Nodes = cell2mat(cellfun(@str2num,Nodes,'UniformOutput',false));
% Nodes = dlmread(fname,',',[idx1,0,idx2-2,2]); % another way to do it

% Elements: number 1st 2nd 3rd 4th (number of elements is where is data)
Elements = AllText(idx2+1:idx3(1)-1) ;
Elements = cell2mat(cellfun(@str2num,Elements,'UniformOutput',false));
end