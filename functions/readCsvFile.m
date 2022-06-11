%% Ler Arquivo
% 
%  Autor: Marcos Netto
%  Email: mpnetto88@gmail.com
%% Retorna Conteudo de arquivo
% Ler conteudo dentro de um arquivo e retorna em formato de matrix

function [Data,Labels] = readCsvFile(file)
    % Abre Arquivo
    fac = fopen(file);
    
    % Pega o conteudo do cabeçalho
    cabecalho = fgets(fac);
    
    % Retorna ponteiro a posicao inicial
    frewind(fac);
    
    % Estipula delimitador (Nesse caso o tab)
    %delimiter = char(9);
    delimiter =  ',';
     
    % Pega a quantidade de Colunas
    numCols = numel(strfind(cabecalho,delimiter));
    
    % Checa se o ultimo caracter é um delimitador
    if cabecalho(end-1) ~= delimiter
        numCols = numCols + 1;
    end

    formatLabel = repmat('%s',1,numCols);
    formatData = '%s%s%f%s';

    conteudo = textscan(fac, formatData, 'Delimiter', delimiter, 'HeaderLines', 1, 'ReturnOnError', false);
    
    data = conteudo(:,3);
    conteudo(:,3) = [];
    %Fecha o arquivo
    fclose(fac);

    conteudo = [conteudo{1:end}];
    data = [data{1:end}];
    
    [~,n,i] = unique(conteudo(:,1));
    
    Labels = conteudo(n, 3)';
    

    for k = 1 : max(i)
      Data(:,k) = data(i == k);
    end
end