classdef Project < handle
    properties
        Group         % Group object containing the project's data
        Configuration % Configuration object containing the project's settings
    end

    methods
        % Constructor for Project class.
        function project = Project(configuration)
            project.Configuration = configuration;
            project.Group = Group(configuration);
        end

        % Run the project's calculation process.
        function run(self)
            % Measure the total time taken for calculations
            totalTimeStart = tic;

            % Calculate the graph features
            self.Group.calculate();

            % Calculate the graph measures
            self.Group.calculateMeasures();

            % Write the calculated measures to output files
            self.Group.write();

            % Display the total calculation time
            totalTimeEnd = toc(totalTimeStart);
            disp("Total calculation time is: " + totalTimeEnd + " seconds.");
        end
    end
end