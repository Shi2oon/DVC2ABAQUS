function fileout = pythonFileName(filein)

count=0;
count=count+1;
[folder{count},file{count},ext] = fileparts(filein);
for i=1:100
count=count+1;  
[folder{count},file{count}]= fileparts(folder{count-1});
if strcmp(folder{count} ,folder{count-1}) 
    break;
end
end

fileout = [folder{end} '\' file{count-1}];
for ii=i-1:-1:1
    fileout = [fileout '\\' file{ii}];
end

fileout = [fileout ext];