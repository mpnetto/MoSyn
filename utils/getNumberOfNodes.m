function numberOfNodes = getNumberOfNodes(file)
    % getNumberOfNodes extracts the number of nodes from a text file.
    %
    % @param file: Input file name with path (file should be of text type).
    % @return numberOfNodes: The number of nodes in the file.

    % Detect the import options for the input file.
    opts = detectImportOptions(file, 'FileType', 'text');

    % Retrieve the number of nodes (column headers) from the import options.
    numberOfNodes = length(opts.VariableNames);

end