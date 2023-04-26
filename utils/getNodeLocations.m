function XY = getNodeLocations(labels, filename)
    % getNodeLocations extracts the node locations from a given file.
    %
    % @param labels: Array of node labels.
    % @param filename: Input file name with path.
    % @return XY: A matrix containing the X and Y coordinates of each node.

    delimiter = '\t';
    startRow = 2;
    endRow = inf;

    formatSpec = '%s%f%f%[^\n\r]';

    fileID = fopen(filename, 'r');

    dataArray = textscan(fileID, formatSpec, endRow(1) - startRow(1) + 1, 'Delimiter', delimiter, 'HeaderLines', startRow(1) - 1, 'ReturnOnError', false);

    fclose(fileID);

    X = zeros(1, length(labels));
    Y = zeros(1, length(labels));

    for i = 1:length(labels)
        id = find(strcmp(upper(dataArray{1, 1}), upper(labels(i))), 1, 'first');
        
        if isempty(id)
            X(i) = 0.0;
            Y(i) = 0.0;
        else
            X(i) = dataArray{1, 2}(id);
            Y(i) = dataArray{1, 3}(id);
        end
    end

    XY = [X' Y'];
end