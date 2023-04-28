% MeasureRegistry is a singleton class responsible for managing the
% registration and retrieval of graph measures. It automatically discovers
% and registers Measure classes found in the 'measures' directory. The
% class provides methods to perform calculations on time-varying graphs and
% access the calculated data.
%
% Usage:
%   registry = MeasureRegistry.getInstance();
%   measures = registry.getMeasures();

classdef MeasureRegistry

    properties (Access = private)
        measures % Container for registered measures
    end

    methods (Access = private)
        % Constructor for MeasureRegistry
        function obj = MeasureRegistry()
            obj.measures = containers.Map();
            measureFiles = dir(strcat([fileparts(which('mosyn')) filesep 'measures'], '/*.m'));
            for i = 1:length(measureFiles)
                obj.registerMeasureFromFile(fullfile(measureFiles(i).folder, measureFiles(i).name));
            end
        end

        % Register a measure from a file.
        %
        % @param filename: The filename of the measure class file.
        function registerMeasureFromFile(obj, filename)
            [~, name, ~] = fileparts(filename);
            className = str2func(name);
            measure = className();
            obj.registerMeasure(measure);
        end
    end

    methods (Static)
        % Define a static method to get the instance of the class.
        function obj = getInstance()
            persistent uniqueInstance

            if isempty(uniqueInstance)
                uniqueInstance = MeasureRegistry();
            end

            obj = uniqueInstance;
        end
    end

    methods
        % Register a measure object in the measures container.
        %
        % @param measure: The Measure object to register.
        function registerMeasure(obj, measure)
            obj.measures(measure.name) = measure;
        end

        % Retrieve a measure object by its name.
        %
        % @param name: The name of the measure to retrieve.
        % @return measure: The Measure object with the specified name.
        function measure = getMeasure(obj, name)
            measure = obj.measures(name);
        end

        % Perform calculations for all registered measures on the given time-varying graph.
        %
        % @param tvg: The TimeVaryingGraph object to perform calculations on.
        function calculateMeasures(obj, tvg)
            keys = obj.measures.keys;

            for i = 1:length(keys)
                measure = obj.measures(keys{i});
                measure.calculateMeasure(tvg);
            end
        end

        % Get an array of all registered Measure objects.
        %
        % @return measuresArray: An array containing all registered Measure objects.
        function measuresArray = getMeasures(obj)
            measuresArray = Measure.empty();

            measures = obj.measures.values;
            for i = 1:length(measures)
                measuresArray(i) = measures{i};
            end
        end

        % Get properties of all registered measures as a cell array.
        %
        % @param propertiesNames: A cell array of property names to retrieve.
        % @return measuresProperties: A cell array containing the specified properties for each measure.
        function measuresProperties = getMeasuresProperties(obj, propertiesNames)
            measures = obj.getMeasures;

            for i = 1:length(propertiesNames)
                measuresProperties(:,i) =  {measures.(propertiesNames{i})}';
            end
        end

        % Get an array of active Measure objects that support the specified graph type.
        %
        % @param graphType: The graph type to filter active measures by.
        % @return measuresArray: An array containing
         % active Measure objects that support the specified graph type.
        function measuresArray = getActiveMeasures(obj, graphType)
            measuresArray = Measure.empty();
    
            measures = obj.measures.values;
    
            for i = 1:length(measures)
                measure = measures{i};
    
                if (measure.active && ismember(graphType, measure.supportedGraphs))
                    measuresArray(end+1) = copy(measure);
                end
            end
        end

        % Get the names of all active measures.
        %
        % @return measureNames: A cell array containing the names of all active measures.
        function measureNames = getMeasuresNames(obj)
            measureNames = {obj.getActiveMeasures.name};
        end
    
        % Write the calculated data for all active measures.
        function writeData(obj)
            for j = 1:length(obj.getActiveMeasures)
    
                measure = obj.getActiveMeasures{j};
    
                measure.write()
            end
        end
    end
end