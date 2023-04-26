function writeFile(matrix, head, file, ext)
    % writeFile converts a given matrix into a table and writes it to a file.
    %
    % @param matrix: Input data in cell or array format.
    % @param head: Cell array of header (column) names for the table.
    % @param file: File name (without extension) where the table will be written.
    % @param ext: File extension (including the dot) for the output file.

    % Convert input matrix to a table.
    if iscell(matrix)
        table = cell2table(matrix);
    else
        table = array2table(matrix);
    end

    % Construct the output file name.
    fileName = [file ext];

    % Set table headers.
    table.Properties.VariableNames = head;

    % Write the table to the file with tab delimiter.
    writetable(table, fileName, 'Delimiter', '\t');
end