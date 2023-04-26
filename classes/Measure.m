% Measure is an abstract class serving as a base for creating specific
% graph measures. It defines a common structure and properties for all
% derived classes. The main purpose of this class is to calculate a
% specific graph measure and perform basic data manipulation, such as
% computing averages, standard deviations, and variation coefficients.
%
% Usage:
%   Derived classes should inherit from this class and implement the
%   abstract calculate methods (calculateBinary and calculateWeighted).

classdef (Abstract) Measure < handle & matlab.mixin.Copyable & matlab.mixin.Heterogeneous
    properties (Abstract)
        active          % (logical) Active status of the measure
        name            % (string) Name of the measure
        calculateType   % (string) Type of calculation (e.g., 'Frame')
        supportedGraphs % (cell array of strings) List of supported graph types
    end

    properties
        Graph                % (Graph object) Graph for which the measure is calculated
        value                % (numeric array) Measure values
        frameAverage         % (numeric) Average of the values per frame
        average              % (numeric) Overall average of the values
        standardDeviation    % (numeric) Standard deviation of the values
        variationCoefficient % (numeric) Variation coefficient of the values
    end

    methods (Abstract)
        % Abstract method for calculating the binary measure for a given
        % adjacency matrix.
        %
        % @param adjMatrix: (numeric matrix) A square adjacency matrix representing the
        %                   graph of interest.
        calculateBinary(adjMatrix);

        % Abstract method for calculating the weighted measure for a given
        % adjacency matrix.
        %
        % @param adjMatrix: (numeric matrix) A square adjacency matrix representing the
        %                   graph of interest.
        calculateWeighted(adjMatrix);
    end

    methods
        % Calculate the measure based on the graph type (binary or weighted)
        % and graph dimensions.
        %
        % @param graph: (Graph object) The graph for which the measure is calculated.
        function calculate(self, graph)

            self.Graph = graph;

            switch graph.Type
                case 'GraphBD'
                    calculateFun = @self.calculateBinary;
                case 'GraphWD'
                    calculateFun = @self.calculateWeighted;
                otherwise
                    error('Unsupported graph type');
            end

            if(self.Graph.Dimensions == 3)
                self.calculateTensor(graph.Data, calculateFun)
            else
                self.calculateGraph(graph.Data, calculateFun)
            end
        end

        % Calculate the measure for a 2D graph and save the result.
        %
        % @param graphData: (numeric matrix) The 2D graph data.
        % @param calculateFun: (function handle) The measure calculation function.
        function calculateGraph(self, graphData, calculateFun)
            values = calculateFun(graphData);

            self.saveAllData(values)
        end

        % Calculate the measure for a 3D graph (tensor) and save the result.
        %
        % @param graphData: (numeric tensor) The 3D graph data.
        % @param calculateFun: (function handle) The measure calculation function.
        
        function calculateTensor(self, graphData, calculateFun)

            if strcmp(self.calculateType, 'Frame')

                parfor i = 1:length(graphData)
    
                    AdjacencyMatrix = graphData(:,:,i);
    
                    % Compute the measure values for the current graph.
                    values(:,i) = calculateFun(AdjacencyMatrix);

                end

            else
                
                 values = calculateFun(graphData);
            end
            % Save the calculated measure values for all graphs.
            self.saveAllData(values)
        end

        % Calculate the graph measure for a given adjacency matrix and
        % save the result.
        %
        % @param adjMatrix: A square adjacency matrix representing the
        %                   graph of interest.
        function calculateMeasure(self, adjMatrix)
            val = self.calculate(adjMatrix);
            self.saveData(val)
        end

        % Save the calculated measure value to the object's property.
        %
        % @param value: The calculated measure value.
        function saveData(self, value)
            self.value = [self.value ; value];
        end

        % Save all data for the measure, including frame average, overall
        % average, standard deviation, and variation coefficient.
        %
        % @param value: The calculated measure value.
        function saveAllData(self, value)
            self.value = value;

            if(~isempty(self.Graph) && self.Graph.Dimensions == 3)
                self.frameAverage = mean(value)';
                self.average = mean(self.frameAverage);
                self.standardDeviation = std(self.frameAverage ,1);
                self.variationCoefficient = self.standardDeviation / self.average;
            else
                self.average = mean(self.value);
                self.standardDeviation = std(self.value ,1);
                self.variationCoefficient = self.standardDeviation / self.average;
            end

           
        end

       

    end
    methods (Sealed)
        % Get the value of the specified property.
        %
        % @param propertyName: (string) The name of the property.
        % @return value: (varies depending on the property) The value of the specified property.
        function value = getPropertyValue(self, propertyName)
            value = eval(strcat('self.', propertyName));
        end

        % Get the data for the specified properties.
        %
        % @param propertiesNames: (cell array of strings) The names of the properties.
        % @return measureData: (cell array) The data for the specified properties.
        function measureData = getMeasureData(self, propertiesNames)

            for i = 1:length(propertiesNames)
                measureData(i) = {self.getPropertyValue(propertiesNames{i})};
            end

        end

        % Write the frame average of the measure to a file.
        %
        % @param path: (string) The path where the file should be saved.
        % @param name: (string) The name of the file (without extension).
        function write(self, path, name)

            len = length(self.frameAverage)-1;

            writeFile([(0:1:len)', self.frameAverage], {'Frame', self.name}, [path '\' name '_' self.name], '.txt')
        end
    end
end