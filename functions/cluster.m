%% Weighted Degree
% 
%  Autor: Marcos Netto
%  Email: mpnetto88@gmail.com
%% Calculate 
% implementar comentario

function [cl, clNode] = cluster(A)
    
    numNodes = length(A);
    
    clNode=zeros(numNodes,1);

    for node=1:numNodes
        neigh=find(A(node,:));
        KK=length(neigh);
        if KK>=2
            S=A(neigh,neigh);
            clNode(node)=sum(S(:))/(KK*(KK-1));
        end
    end
    
    cl = mean(clNode);
end
