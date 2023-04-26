% PathLength class computes the average shortest path length obetween all 
% pairs of nodes in a given adjacency matrix. This class inherits from the 
% Measure class.
%
% The shortest path length is the minimum number of edges that must be 
% traversed to get from one node to another in a graph.
%
% Usage:
%   pl = PathLength();
%   values = pl.calculate(adjMatrix);

classdef PathLength < Measure
    properties
        name = 'PathLength';            % (string) Name of the measure
        active = false;                 % (logical) Active status of the measure
        calculateType = 'Frame';        % (string) Type of calculation (e.g., 'Frame', 'All')
        supportedGraphs = ["GraphBD"];  % (cell array of strings) List of supported graph types
    end

    methods
        
         % Calculate the average shortest path length for the binary graph
        % represented by the given adjacency matrix.
        %
        % @param adjMatrix: The adjacency matrix of the graph.
        % @return value: A column vector containing the average shortest
        %                path length for each node.
        function value = calculateBinary(~, adjMatrix)
          % Determine the number of nodes in the graph.
            n = length(adjMatrix);

            % Initialize the distance matrix with infinite values.
            d = inf(n);

            % Fill in the distance matrix with edge weights for adjacent nodes.
            for i = 1:n
                for j = 1:n
                    if adjMatrix(i, j) ~= 0 % If there is an edge between i and j
                        d(i, j) = double(adjMatrix(i, j));
                    end
                end
            end

            % Implement the Floyd-Warshall algorithm to find the shortest path
            % length between all pairs of nodes.
            for k = 1:n
                for i = 1:n
                    for j = 1:n
                        % Update the distance matrix if a shorter path is found.
                        if d(i, j) > d(i, k) + d(k, j)
                            d(i, j) = d(i, k) + d(k, j);
                        end
                    end
                end
            end

            % Remove the diagonal elements (self-distances) from the distance matrix.
            d(1:1+length(d):end) = 0;

            % Initialize the output column vector.
            value = zeros(n, 1);

            % Compute the average shortest path length for each node.
            for i = 1:n
                totalPaths = sum(d(i, ~isinf(d(i, :))));
                numPaths = length(nonzeros(d(i,:)~=Inf)); % Exclude infinite distances and diagonal
                value(i) = totalPaths / numPaths;
            end
        end

        % Calculate the average shortest path length for the weighted graph
        % represented by the given adjacency matrix.
        %
        % @param adjMatrix: The adjacency matrix of the graph.
        % @return value: A column vector containing the average shortest
        %                path length for each node.
        function value = calculateWeighted(~, adjMatrix)
            value = calculatePathLength(adjMatrix);
        end

        % Calculate the average shortest path length for the graph
        % represented by the given adjacency matrix (common function for
        % both binary and weighted graphs).
        %
        % @param adjMatrix: The adjacency matrix of the graph.
        % @return value: A column vector containing the average shortest
        %                path length for each node.
        function value = calculatePathLength(adjMatrix)

            
        end
    end
end