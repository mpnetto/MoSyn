% ClusterCoefficient class computes the clustering coefficient of each node
% in a given adjacency matrix. This class inherits from the Measure class.
%
% The clustering coefficient is a measure of the degree to which nodes
% in a graph tend to cluster together.
%
% Usage:
%   cc = ClusterCoefficient();
%   values = cc.calculate(adjMatrix);

classdef ClusterCoefficient < Measure
    properties
        name = 'ClusterCoefficient';         % (string) Name of the measure
        active = true;                       % (logical) Active status of the measure
        calculateType = 'Frame';             % (string) Type of calculation (e.g., 'Frame', 'All')
        supportedGraphs = ["GraphBD", "GraphWD"]; % (cell array of strings) List of supported graph types
    end
    
    methods
        % Calculate the binary clustering coefficient for each node in the 
        % given adjacency matrix.
        %
        % @param adjMatrix: (numeric matrix) The adjacency matrix of the graph.
        % @return value: (numeric column vector) The clustering coefficient for each node.
        function value = calculateBinary(~, adjMatrix)
           % Determine the number of nodes in the graph.
            numNodes = length(adjMatrix);
            
            % Initialize the output column vector with zeros.
            value = zeros(numNodes,1);
            
            % Loop through each node in the graph.
            for node = 1:numNodes
                % Find the neighbors of the current node.
                neigh = find(adjMatrix(node,:));
                
                % Compute the number of neighbors (degree) for the current node.
                kk = length(neigh);
                
                % Calculate the clustering coefficient for the current node
                % if it has at least two neighbors.
                if kk >= 2
                    % Extract the subgraph induced by the neighbors.
                    subGraph = adjMatrix(neigh,neigh);
                
                    % Compute the clustering coefficient for the current node.
                    value(node) = sum(subGraph(:))/(kk*(kk-1));
                end
            end
        end

        % Calculate the weighted clustering coefficient for each node in the 
        % given adjacency matrix.
        %
        % @param adjMatrix: (numeric matrix) The adjacency matrix of the graph.
        % @return value: (numeric column vector) The weighted clustering coefficient for each node.
        function value = calculateWeighted(~, adjMatrix)
             % Determine the number of nodes in the graph.
            numNodes = length(adjMatrix);

            % Initialize the output column vector with zeros.
            value = zeros(numNodes,1);

            % Loop through each node in the graph.
            for node = 1:numNodes
                % Find the neighbors of the current node.
                neigh = find(adjMatrix(node,:));

                % Compute the number of neighbors (degree) for the current node.
                kk = length(neigh);

                % Calculate the weighted clustering coefficient for the current node
                % if it has at least two neighbors.
                if kk >= 2
                    % Extract the subgraph induced by the neighbors.
                    subGraph = adjMatrix(neigh, neigh);

                    % Calculate the strength of the current node.
                    strength = sum(adjMatrix(node, :));

                    % Calculate the weighted clustering coefficient for the current node.
                    T_i = sum(subGraph(:));
                    value(node) = (2 * T_i) / (kk * (kk - 1) * strength);
                end
            end
        end
    end
end
