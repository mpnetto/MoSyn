classdef TVG < GraphFeatures

    methods

        % Constructor for TVG class.
        function obj = TVG(binaryTVG)
            obj.Data = binaryTVG;
            obj.Type = 'GraphBD';
            obj.Dimensions = 3;

            % Obtain the singleton instance of the MeasureRegistry class.
            measureRegistry = MeasureRegistry.getInstance();

            % Retrieve the list of measures from the measure registry.
            obj.Measures = measureRegistry.getActiveMeasures(obj.Type);
        end

        % Calculate various graph measures for each binary time-varying graph
        % and saves the calculated values.
        function calculateMeasures(obj)
            multiWaitbar('Calculation measures', 0, 'Color', 'g', 'CanCancel', 'on', 'CancelFcn', @(a, b) disp(['Cancel ', a]));

            % Loop through each measure in the measure registry.
            for i = 1:length(obj.Measures)
                tic
                % Extract the current measure.
                measure = obj.Measures(i);
                disp(measure.name);

                % Calculate the measure for the current TVG object.
                measure.calculate(obj);
                abort = multiWaitbar('Calculation measures', i / length(obj.Measures));
                toc
            end

            multiWaitbar('Calculation measures', 'Close');
        end

        % Get the graph at a specific position in the time-varying graph data.
        function gr = getGraphAtPosition(obj, pos)
            gr = graph(obj.Data(:, :, pos));
        end

        % Write the calculated measures to output files.
        function write(self)
            self.Measures.write();
        end
    end
end
