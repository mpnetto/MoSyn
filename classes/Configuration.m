% Configuration is a class that holds the configuration settings for a
% given dataset. It contains various properties related to the dataset,
% such as the path, files, number of nodes, node location, and parameters
% for data processing.

classdef Configuration < handle & matlab.mixin.Copyable
    properties
       Path                % The path to the dataset directory
       Files               % The files within the dataset directory
       Filename            % The name of the dataset file
       Name                % The name of the dataset
       Extension           % The file extension of the dataset
       Size                % The size of the dataset
       NodesLabels         % The labels of the nodes in the dataset
       NumberOfNodes       % The number of nodes in the dataset
       NodeLocation        % The location of the nodes in the dataset
       InitialTime         % The initial time of the data processing window
       FinalTime           % The final time of the data processing window
       RealTime            % Indicates if the real-time data is being used
       InitialFrequency    % The initial frequency for processing
       FinalFrequency      % The final frequency for processing
       HasFrequency        % Indicates if the dataset has frequency data
       RemoveTimeColumn    % Indicates if the time column should be removed
       AcqRate             % The acquisition rate of the dataset
       SlidWindow          % The sliding window size for processing
       TaoMin              % The minimum Tao value for processing
       TaoMax              % The maximum Tao value for processing
       Threshold           % The threshold value for processing
       TscLim              % The TSC limit value for processing
    end

    properties (Dependent)
        RealEndTime        % The real end time of the data processing window
    end

    methods
        % Constructor for the Configuration class
        function configuration = Configuration(path, files, nodeLocation, initialTime, finalTime, realTime, initialFrequency, finalFrequency, removeTimeColumn, acqRate, slidWindow, taoMin, taoMax, threshold, tscLim)
            configuration.Path = path;
            configuration.Files = files;
            configuration.NodeLocation = nodeLocation;
            configuration.InitialTime = initialTime;
            configuration.FinalTime = finalTime;
            configuration.RealTime = realTime;
            configuration.InitialFrequency = initialFrequency;
            configuration.FinalFrequency = finalFrequency;
            configuration.RemoveTimeColumn = removeTimeColumn;
            configuration.AcqRate = acqRate;
            configuration.SlidWindow = slidWindow;
            configuration.TaoMin = taoMin;
            configuration.TaoMax = taoMax;
            configuration.Threshold = threshold;
            configuration.TscLim = tscLim;
        end

        % Get the real end time of the data processing window
        function endTime = get.RealEndTime(obj)
            realEndTime = obj.Size - obj.SlidWindow - obj.TaoMax - 1;

            if (obj.RealTime)
                endTime = realEndTime;
            elseif (obj.FinalTime > realEndTime)
                endTime = realEndTime;
                msg = ['Tamanho do arquivo menor do que o tempo final informado: Usando o tamanho como tempo final. Pressione OK para continuar. Arquivo:' obj.Filename];
                uiwait(msgbox(msg, 'Warning', 'modal'));
            else
                endTime = obj.FinalTime;
            end
        end
    end
end