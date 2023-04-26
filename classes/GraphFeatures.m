% GraphFeatures is an abstract class that represents the features of a
% graph, including its data, type, dimensions, and measures. This class
% provides methods to access the calculated data and measures associated
% with the graph.
%
% Subclasses should implement the 'write' method to handle the output
% of the calculated data.

classdef (Abstract) GraphFeatures < handle
    properties
        Data       % The graph data
        Type       % The type of the graph (e.g., 'GraphBD', 'GraphWD')
        Dimensions % The dimensions of the graph
        Measures   % A list of Measure objects associated with the graph
    end

    methods (Abstract)
        % Write the calculated data for the graph features.
        write(self);
    end

    methods
        % Get the data of the specified properties for all measures.
        %
        % @param propertiesNames: A cell array containing the names of the
        %                         properties to retrieve.
        % @return measuresData: A cell array containing the requested data
        %                      for all measures.
        function measuresData = getMeasuresData(obj, propertiesNames)
            measuresData = {};

            for i = 1:length(obj.Measures)
                measuresData(end+1,:) = obj.Measures(i).getMeasureData(propertiesNames);
            end
        end

        % Get a specific measure by its name.
        %
        % @param measureName: The name of the measure to retrieve.
        % @return measure: The Measure object with the specified name, or
        %                 an empty array if not found.
        function measure = getMeasure(obj, measureName)
            for i = 1:length(obj.Measures)
                if strcmp(obj.Measures{i}.name, measureName)
                    measure = obj.Measures{i};
                    return;
                end
            end

            measure = [];
        end
    end
end