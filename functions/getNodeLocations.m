function  XY = getNodeLocations(labels, filename) 
%% Initialize variables.

delimiter = '\t';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format string for each line of text:
formatSpec = '%s%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

X = zeros(1, length (labels));
Y = zeros(1, length (labels));

%% Create output variable
 for i=1:length (labels)
     id=find(strcmp(upper(dataArray{1,1}),upper(labels(i))),1,'first');
     if isempty(id)
         X(i)=0.0;Y(i)=0.0;
     else
        X(i)=dataArray{1,2}(id);
        Y(i)=dataArray{1,3}(id);
     end
 end
 XY=[X' Y'];
end