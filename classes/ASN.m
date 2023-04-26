classdef ASN < GraphFeatures
    properties
    end

    methods

        % Constructor for ASN class.
        function obj = ASN(binaryTVG)
            obj.Data = sum(binaryTVG, 3);
            obj.Type = 'GraphWD';
            obj.Dimensions = 2;

            % Get instance of MeasureRegistry and retrieve active measures.
            measureRegistry = MeasureRegistry.getInstance();
            obj.Measures = measureRegistry.getActiveMeasures(obj.Type);
        end

        % Calculate measures for the ASN object.
        function calculateMeasures(obj)
            multiWaitbar('Calculating ASN ...', 0, 'Color', 'g', 'CanCancel', 'on', 'CancelFcn', @(a, b) disp(['Cancel ', a]));

            % Iterate through the measures and calculate them.
            for i = 1:length(obj.Measures)
                measure = obj.Measures(i);
                measure.calculate(obj);
                abort = multiWaitbar('Calculating ASN ...', i / length(obj.Measures));
            end

            multiWaitbar('Calculating ASN ...', 'Close');
        end

        % Write results to output files.
        function write(self, configuration)
            self.writeAsnCorrelationWithoutZeros(configuration);
            self.writeAsnCorrelationWithoutZeros(configuration);
        end

        % Write ASN correlation without zeros to the output file.
        function writeAsnCorrelationWithoutZeros(self, configuration)
            asnIn = tril(self.Data, -1);

            [target, source, weight] = find(asnIn);
            target = configuration.NodesLabels(target);
            source = configuration.NodesLabels(source);

            type = cell(length(target), 1);
            type(:) = {'Undirected'};

            table = [source', target', num2cell(weight), type];

            writeFile(table, {'Source', 'Target', 'Weight', 'Type'}, [configuration.Path '\' configuration.Name '_ASN_G'], '.txt');
        end

        % Write ASN correlation to the output file.
        function writeAsnCorrelation(self, path, name, nodesLabels)
            asnOut = triu(self.Data + 1, 1);

            [target, source, weight] = find(asnOut);
            weight = weight + 1;

            target = configuration.NodesLabels(target);
            source = configuration.NodesLabels(source);

            type = cell(length(target), 1);
            type(:) = {'Undirected'};

            table = [source', target', num2cell(weight), type];

            writeFile(table, {'Source', 'Target', 'Weight', 'Type'}, [configuration.Path '\' configuration.Name '_ASN_T'], '.txt');
        end
    end
end