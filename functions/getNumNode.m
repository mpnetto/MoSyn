%% Ler Arquivo
% 
%  Autor: Marcos Netto
%  Email: mpnetto88@gmail.com
%% Retorna Conteudo de arquivo
% Ler conteudo dentro de um arquivo e retorna em formato de matrix

function numNodes = getNumNode(file)
    % Open File
    fac = fopen(file);
    
    % Get head content
    head = fgets(fac);
    
    % Set delimiter
    delimiter = char(9);
    
    % Pega a quantidade de Colunas
    numNodes = numel(strfind(head,delimiter));

    % Checa se o ultimo caracter é um delimitador
    if head(end-1) ~= delimiter
        numNodes = numNodes + 1;
    end
    
    %Close File
    fclose(fac);

end