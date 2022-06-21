%% Ler Arquivo
% 
%  Autor: Marcos Netto
%  Email: mpnetto88@gmail.com
%% Retorna Conteudo de arquivo
% Ler conteudo dentro de um arquivo e retorna em formato de matrix

function Data = readLocation(tvg)
    
    if isempty(tvg.NodesLocation)
        for i=1:length(tvg.NodesLabels)
            x(i)=rand(); y(i)=rand();
        end
    else
        lb=readtable(tvg.NodesLocation);
        for i=1:length(tvg.NodesLabels)
            id=find(ismember(upper(lb{:,1}),upper(tvg.Labels{i})));
            if isempty(id)
                x(i)=0;y(i)=0;
            else
                x(i)=lb{id,2};
                y(i)=lb{id,3};
            end
        end
    end
    
    
    Data = [x' y'];
  
end