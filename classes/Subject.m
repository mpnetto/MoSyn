% Subject is a class representing a subject in a study. It holds the
% subject's data, configuration, time-varying graph (TVG), and added
% static network (ASN). The class provides methods to load data, perform
% calculations on the data, and calculate graph measures.
%
% Usage:
%   subject = Subject(configuration, file);

classdef Subject < handle

    properties
       Configuration % Configuration object for the subject
       TVG           % Time-varying graph (TVG) of the subject
       ASN           % Adjacency sequence network (ASN) of the subject
    end

    properties (Dependent)
        Data % Dependent property for subject data
    end

    methods
        % Getter for the Data dependent property.
        function data = get.Data(obj)

            % Detect import options for the file specified in the configuration
            opts = detectImportOptions(obj.Configuration.Filename, 'FileType','text');

            % Store the variable names from the import options as labels
            obj.Configuration.NodesLabels = opts.VariableNames;
            
            if verLessThan('matlab', '9.6')
                data = readFile(obj.Configuration.Filename, obj.Configuration.Extension);
            else
                  % Read the matrix from the file using the import options
                data = readmatrix(obj.Configuration.Filename, opts);
             
            end


            % Remove the first column if 'Remove Time Column' is checked
            if(obj.Configuration.RemoveTimeColumn == "on")
                data(:,1) = [];
                obj.Configuration.NodesLabels(1,:) = [];
            end

            % Apply bandpass filter if 'Allow Band Pass' is checked
            if(obj.Configuration.HasFrequency)
                data = bandpass(data,[obj.Configuration.InitialFrequency obj.Configuration.FinalFrequency], obj.Configuration.AcqRate);
            end

            % Set the size and number of nodes in the configuration
            obj.Configuration.Size = length(data);
            obj.Configuration.NumberOfNodes = length(opts.VariableNames);
        end
    end

    methods
        % Constructor for the Subject class.
        %
        % @param configuration: Configuration object for the subject.
        % @param file: Filename of the subject data file.
        function subject = Subject(configuration, file)
            if nargin > 0
                filename = fullfile(configuration.Path, file);
    
                [~, name, extension] = fileparts(filename);

                subject.Configuration = configuration;
                subject.Configuration.Filename = filename;
                subject.Configuration.Name = name;
                subject.Configuration.Extension = extension;
            end
        end

        % Perform calculations on the subject data to generate TVG and ASN.
        function calculate(obj)

            data = obj.Data;
            initialTime = 0;
            realEndTime = obj.Configuration.RealEndTime;
            binaryTVG = motif_sync_mex(data, initialTime, realEndTime, obj.Configuration.NumberOfNodes, obj.Configuration.SlidWindow, obj.Configuration.TaoMin, obj.Configuration.TaoMax, obj.Configuration.Threshold);
            
            obj.TVG = TVG(binaryTVG);
            obj.ASN = ASN(binaryTVG);
        end

        % Calculate measures for TVG and ASN.
        function calculateMeasures(obj)
            obj.TVG.calculateMeasures()
            obj.ASN.calculateMeasures()
        end

        function write(self)
            path = self.Configuration.Path;
            name = self.Configuration.Name;

            measures = [self.TVG.Measures];

            for i=1:length(measures)
                measures(i).write(path, name)
            end

            self.ASN.write(self.Configuration)

        end

    end

    methods (Static)
        % Create a zeros-initialized Subject array.
        %
        % @param varargin: Dimensions of the Subject array.
        % @return z: A zeros-initialized Subject array.
        function z = zeros(varargin)
            if(nargin == 0)
                z = Subject;
            elseif any([varargin{:}] <= 0)
                z = Subject.empty(varargin{:});
            else
                z = repmat(Subject,varargin{:});
            end
        end
    end

end