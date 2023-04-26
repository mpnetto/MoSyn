classdef Group < handle
    properties
        Subjects % Array of Subject objects in the group
        Measures % Array of Measure objects for the group
    end

    methods
        % Constructor for the Group class.
        function group = Group(configuration)
            measureRegistry = MeasureRegistry.getInstance();
            group.Measures = measureRegistry.getActiveMeasures("GraphBD");
            group.Subjects = zeros(1, length(configuration.Files), 'Subject');

            for i = 1:length(configuration.Files)
                file = configuration.Files{i};
                group.Subjects(i) = Subject(copy(configuration), file);
            end
        end

        % Calculate graph features for all subjects in the group.
        function calculate(obj)
            multiWaitbar('Calculating TVG', 0, 'Color', 'g', 'CanCancel', 'on', 'CancelFcn', @(a, b) disp(['Cancel ', a]));

            for i = 1:length(obj.Subjects)
                obj.Subjects(i).calculate();
                abort = multiWaitbar('Calculating TVG', i / length(obj.Subjects));
            end
            multiWaitbar('CloseAll');
        end

        % Calculate graph measures for all subjects in the group.
        function calculateMeasures(obj)
            multiWaitbar('Subjects', 0, 'Color', 'g', 'CanCancel', 'on', 'CancelFcn', @(a, b) disp(['Cancel ', a]));

            for i = 1:length(obj.Subjects)
                obj.Subjects(i).calculateMeasures();
                abort = multiWaitbar('Subjects', i / length(obj.Subjects));
            end

            multiWaitbar('CloseAll');

            tvgs = [obj.Subjects.TVG];
            measures = [tvgs.Measures];

            for i = 1:length(obj.Measures)
                measure = obj.Measures(i);
                measuresNames = {measures.name};
                logicalMeasures = strcmp(measuresNames, measure.name);
                value = [measures(logicalMeasures).average];
                measure.saveAllData(value);
            end
        end

        % Get measures data for the specified properties.
        % @param propertiesNames: Cell array of property names.
        % @return measuresData: Cell array of measures data.
        function measuresData = getMeasuresData(obj, propertiesNames)
            measuresData = {};

            for i = 1:length(obj.Measures)
                measuresData(end + 1, :) = obj.Measures(i).getMeasureData(propertiesNames);
            end
        end

        % Write the calculated measures to output files.
        function write(self)
            multiWaitbar('Writing Measures', 0, 'Color', 'b', 'CanCancel', 'on', 'CancelFcn', @(a, b) disp(['Cancel ', a]));

            for i = 1:length(self.Subjects)
                self.Subjects(i).write();
                abort = multiWaitbar('Writing Measures', i / length(self.Subjects));
            end

            multiWaitbar('CloseAll');
        end
    end
end