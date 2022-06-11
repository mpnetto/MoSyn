%% Weighted Degree
% 
%  Autor: Marcos Netto
%  Email: mpnetto88@gmail.com
%% Calculate 
% implementar comentario

function pl = pathLength(A)
    
    dist = distance(A);
    numNodes = length(A);

    cin = zeros(1,numNodes);
    for node = 1:1:numNodes
        distNode = dist(:,node);
        cin(node) = sum(distNode(distNode~=Inf))/length(nonzeros(distNode~=Inf));
    end

    cout = zeros(1,numNodes);
    for node = 1:1:numNodes
        distNode = dist(node,:);
        cout(node) = sum(distNode(distNode~=Inf))/length(nonzeros(distNode~=Inf));
    end

    pl = mean([cin; cout],1);
end
