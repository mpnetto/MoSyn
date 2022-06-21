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
    
    pos = find(isletter(head), 1, 'last');

    head = extractBefore(head,pos);
    
    % Estipula delimitador (Nesse caso o tab)
    if any(head == 9)
        delimiter = char(9);
    elseif any(head == 32)
        delimiter = char(32);
    end
    
    % Pega a quantidade de Colunas
    numNodes = numel(strfind(head,delimiter)) + 1;
    
    %Close File
    fclose(fac);

end