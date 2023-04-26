% The Weight class computes the sum of weights for each node in a given
% weighted adjacency matrix. This class inherits from the Measure class.
%
% The sum of weights for a node is the sum of the weights of all edges
% connected to the node in a weighted graph.
%
% Usage:
% wei = Weight();
% values = wei.calculate(adjMatrix);

classdef Weight < Measure
    properties

        name = 'Weight';                % (string) Name of the measure
        active = true;                  % (logical) Active status of the measure
        calculateType = 'All';          % (string) Type of calculation (e.g., 'Frame', 'All')
        supportedGraphs = ["GraphWD"];  % (cell array of strings) List of supported graph types
    end

    methods

        % This method is not applicable for this class since it only supports
        % weighted graphs. It is left empty and won't be used.
        function calculateBinary(~, ~)
        end

        % Calculate the sum of weights for each node in the given adjacency matrix.
        %
        % @param adjMatrix: The weighted adjacency matrix of the graph.
        % @return value: A column vector containing the sum of weights for each node.
        function value = calculateWeighted(~, adjMatrix)
            value = sum(adjMatrix, 2);
        end
    end
end
