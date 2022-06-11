%% Escreve Arquivo
% 
%  Autor: Marcos Netto
%  Email: mpnetto88@gmail.com
%% Escreve Conteudo em arquivo
% Recebe conteudo em formato de tabela e escreve em arquivo

function writeFile(matrix, head, file, ext)

    if(iscell(matrix))
        table = cell2table(matrix);
    else
        table = array2table(matrix);
    end

     fileName = [file ext];

     table.Properties.VariableNames = head;

     writetable(table,fileName,'Delimiter','\t');
end