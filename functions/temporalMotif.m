
function [g] = temporalMotif(tvg, n, labels)


% motifs de quatro
% gráficos temporais para os 10 que mais se repetem
    numNodes = size(tvg,1);

    nodeMatrix = zeros([numNodes numNodes numNodes]);
    
    logicalArray = arrayfun(@(x) ~~sum(tvg(:,:, x:(n-1)+x),3),1:n:size(tvg,3)-n,'UniformOutput',false);

    logicalMatrix = cat(3,logicalArray{:});

    serie = 1:numNodes;


    for m = 1: size(logicalMatrix,3)
        g = graph(logicalMatrix(:,:,m));    
    
        for i = 1:numNodes
            neighs = neighbors(g,i);
        
            for j = 1:length(neighs)-1
        
                nodes = sort([ neighs(j:j+1); i]);
        
                nodeMatrix(nodes(1),nodes(2),nodes(3)) =  nodeMatrix(nodes(1),nodes(2),nodes(3)) + 1;
        
            end
        
        end
    end

    tvgRea = zeros([numNodes numNodes]);

    for i = serie
        for j = serie
             for k = serie
                 value = nodeMatrix(i,j,k);
                if( value > 0)
                    tvgRea(i,j) =  tvgRea(i,j) + value;
                    tvgRea(j,k) =  tvgRea(j,k) + value;
                end

             end
        end
    end

    g = graph(tvgRea ,'upper');

    coord = getNodeLocations(labels,'Location32.txt');
    x = coord(:,1);
    y = coord(:,2);

    weights =  g.Edges.Weight;
    NormWei = weights./max(weights);
    cm = colormap(winter);
    for yy=1:length(NormWei)
        colorID(yy,1) = max(1, sum(NormWei(yy,1) > [0:1/length(cm(:,1)):1]));
    end
    myColor = cm(colorID, :);

    p = plot(g,  'XData',x,'YData',y, 'LineWidth', (NormWei) * 10, 'EdgeColor', myColor);
    labelnode(p, serie, labels)

end

