%% Calculate TSC
% 
%  Autor: Marcos Netto
%  Email: mpnetto88@gmail.com
%% Calculate 
function [tsc n] = tsc(T, threshold)

time = 0;
Edges = T(:,:,1);
Edges = Edges(:);

limit = threshold * length(Edges);
values = 0;
for i = 2:length(T(1,1,:))
    time = time + 1;
    %Transform the matrix into an edge array
    E = T(:,:,i);
    E = E(:);
    
    % Apply an or function to see how many edges have been activated so far
    Edges = Edges | E;
    value  = nnz(Edges);
    
    %If value of edges activate reached the limit, save the time;
    if value >= limit
        
       Edges = zeros(length(Edges),1);
       values = [values;time];
       time = 0;
    end
end

% remove first value which is -
values(1) = [];


% Return the synchronization times and how many time it reached the limit
tsc = values;
n = length(values);

end

