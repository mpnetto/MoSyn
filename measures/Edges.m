% Edges class computes the number of edges in a given adjacency matrix.
% This class inherits from the Measure class.
%
% The number of edges is the total count of connections between the nodes
% in a graph.
%
% Usage:
%   edges = Edges();
%   value = edges.calculate(adjMatrix);

classdef Edges < Measure
    properties
        name = 'Edges';                         % (string) Name of the measure
        active = true;                         % (logical) Active status of the measure
        calculateType = 'Frame';               % (string) Type of calculation (e.g., 'Frame')
        supportedGraphs = ["GraphBD"];         % (cell array of strings) List of supported graph types
    end

    methods
        % Calculate the number of edges in the given binary adjacency matrix.
        %
        % @param adjMatrix: (numeric matrix) The adjacency matrix of the graph.
        % @return value: (numeric) The number of edges in the graph.
        function value = calculateBinary(~, adjMatrix)
            % Sum the adjacency matrix and divide by 2 to account for
            % double counting of edges.
            value = sum(sum(adjMatrix))/2;
        end

        % Calculate the number of edges in the given weighted adjacency matrix.
        %
        % @param adjMatrix: (numeric matrix) The adjacency matrix of the graph.
        % @return value: (numeric) The number of edges in the graph.
        function value = calculateWeighted(~, adjMatrix)
            % Sum the non-zero values in the adjacency matrix and divide
            % by 2 to account for double counting of edges.
            value = sum(sum(adjMatrix ~= 0))/2;
        end
    end
end
