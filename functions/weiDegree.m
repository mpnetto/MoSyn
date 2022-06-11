%% Weighted Degree
% 
%  Autor: Marcos Netto
%  Email: mpnetto88@gmail.com
%% Calculate 
% implementar comentario

function [degree, mean ,inDegree, outDegree] = weiDegree(A)
    
    inDegree = sum(A,1);
    outDegree = sum(A,2);
    
    degree = inDegree + outDegree;
    
    numNodes = length(A);
    
    mean = sum(degree)/numNodes;
end
