
function [data] = readFile(file, extension)

    fac = fopen(file);
    
    cabecalho = fgets(fac);
    
    frewind(fac);
   
    delimiter = char(9);
    numCols = numel(strfind(cabecalho,delimiter));
    
    if cabecalho(end-1) ~= delimiter
        numCols = numCols + 1;
    end
    
    if (extension == ".csv")
        formatData = '%s%s%f%s';
    else
        formatData = repmat('%f',1,numCols);
    end

    content = textscan(fac, formatData, 'Delimiter', delimiter, 'HeaderLines', 1, 'ReturnOnError', false);
    
    fclose(fac);

    data = [content{1:end}];
end