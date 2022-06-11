%% Ler Arquivo
% 
%  Autor: Marcos Netto
%  Email: mpnetto88@gmail.com
%% Retorna Conteudo de arquivo
% Ler conteudo dentro de um arquivo e retorna em formato de matrix

function [Data,Labels] = readFile(file)
    % Abre Arquivo
    fac = fopen(file);
    
    % Pega o conteudo do cabeçalho
    cabecalho = fgets(fac);
    
    % Retorna ponteiro a posicao inicial
    frewind(fac);
    
    % Estipula delimitador (Nesse caso o tab)
    delimiter = char(9);
    
    % Pega a quantidade de Colunas
    numCols = numel(strfind(cabecalho,delimiter));
    
    % Checa se o ultimo caracter é um delimitador
    if cabecalho(end-1) ~= delimiter
        numCols = numCols + 1;
    end

    formatLabel = repmat('%s',1,numCols);
    formatData = repmat('%f',1,numCols);

    % Pega todo o conteudo do arquivo e coloca em um array
    cabecalho = textscan(fac, formatLabel, 1, 'Delimiter', delimiter, 'ReturnOnError', false);

    conteudo = textscan(fac, formatData, 'Delimiter', delimiter, 'HeaderLines', 1, 'ReturnOnError', false);
    
    %Fecha o arquivo
    fclose(fac);

    Labels = [cabecalho{1:end}]';
    Data = [conteudo{1:end}];
end