% Hub class computes a binary vector that indicates if a node is a hub in a 
% given adjacency matrix. This class inherits from the Measure class.
%
% A hub is a node with a degree significantly higher than the average degree
% of the nodes in the graph (mean + 2 standard deviations).
%
% Usage:
%   hub = Hub();
%   values = hub.calculate(adjMatrix);

classdef Hub < Measure
    properties
        name = 'Hub';                  % (string) Name of the measure
        active = true;                 % (logical) Active status of the measure
        calculateType = 'All';         % (string) Type of calculation (e.g., 'Frame', 'All')
        supportedGraphs = ["GraphBD"]; % (cell array of strings) List of supported graph types
    end

    methods
        % Calculate the hub status for each node in the given adjacency
        % matrix.
        %
        % @param adjMatrix: The adjacency matrix of the graph.
        % @return value: A binary row vector indicating if a node is a hub.

        function value = calculateBinary(obj, adjMatrix)
            % Find the Degree measure object in the graph's measures.
            idx = arrayfun(@(x) strcmp(x.name, 'Degree'), obj.Graph.Measures);

            % Get the Degree measure object.
            degreeMeasure = obj.Graph.Measures(idx);

            % Get the degree values for each node.
            degree = degreeMeasure.value;

            % Calculate the hub status for each node by comparing its degree to
            % the average degree plus 2 standard deviations.
            value = (degree > mean(degree) + 2 * std(degree, 1))';
        end

        % This method is intentionally left empty, as the Hub measure is not
        % applicable to weighted graphs.
        function value = calculateWeighted(~, adjMatrix)
        end
    end
end