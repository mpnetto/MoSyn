%% Weighted Degree
% 
%  Autor: Marcos Netto
%  Email: mpnetto88@gmail.com
%% Calculate 
% implementar comentario

function dist = distance(A)
    
    numNodes = length(A);
    len = 1;  
    dist = zeros(numNodes,numNodes);

    Lpath = A;
    Idx = true;
    while any(Idx(:))
        len = len+1;
        Lpath = Lpath.*A;
        Idx = (Lpath~=0)&(dist==0);
        dist(Idx) = len;
    end

    dist(~dist) = inf;
    dist(1:1+numNodes:end) = 0;

end
