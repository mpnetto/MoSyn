% Degree class computes the degree of each node in a given adjacency 
% matrix. This class inherits from the Measure class.
%
% The degree is the number of edges incident to a node in a graph.
%
% Usage:
%   deg = Degree();
%   values = deg.calculate(adjMatrix);

classdef Degree < Measure
    properties
        name = 'Degree';             % (string) Name of the measure
        active = true;               % (logical) Active status of the measure
        calculateType = 'Frame';     % (string) Type of calculation (e.g., 'Frame', 'All')
        supportedGraphs = ["GraphBD", "GraphWD"]; % (cell array of strings) List of supported graph types
    end

    methods

        % Calculate the binary degree for each node in the given adjacency matrix.
        %
        % @param adjMatrix: (numeric matrix) The adjacency matrix of the graph.
        % @return value: (numeric column vector) The degree for each node.
        function value = calculateBinary(~, adjMatrix)
            % Sum the adjacency matrix along the second dimension to
            % obtain the degree of each node.
            value = sum(adjMatrix, 2);
        end

        % Calculate the weighted degree for each node in the given adjacency matrix.
        %
        % @param adjMatrix: (numeric matrix) The adjacency matrix of the graph.
        % @return value: (numeric column vector) The weighted degree for each node.
        function value = calculateWeighted(~, adjMatrix)
            value = sum(adjMatrix, 2);
        end
    end
end