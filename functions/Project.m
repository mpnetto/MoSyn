classdef Project < handle
   
    properties (Constant)
            
        DEGREE = 1;
        DEGREE_NAME = 'Degree';
        DEGREE_TXT = 'The degree of a project is the average of the subjects degree.';
        
        PATHLENGTH = 2;
        PATHLENGTH_NAME = 'PathLength';
        PATHLENGTH_TXT = 'The path length of a project is the average of the subjects  path length.';

        CLUSTERCOEFFICIENT = 3;
        CLUSTERCOEFFICIENT_NAME = 'ClusterCoefficient';
        CLUSTERCOEFFICIENT_TXT = 'The cluster coefficient of a project is the average of the subjects cluster coefficient.';

        RADIUS = 4;
        RADIUS_NAME = 'Radius';
        RADIUS_TXT = 'The radius is the minimum eccentricity.';
        
        DIAMETER = 5;
        DIAMETER_NAME = 'Diameter';
        DIAMETER_TXT = 'The diameter is the maximum eccentricity.'
        
        BETWEENNESS = 6;
        BETWEENNESS_NAME = 'Betweenness';
        BETWEENNESS_TXT = 'Node betweenness centrality of a node is the fraction of all shortest paths in the graph that contain a given node. Nodes with high values of betweenness centrality participate in a large number of shortest paths.';
        
        CLOSENESS = 7;
        CLOSENESS_NAME = 'Closeness';
        CLOSENESS_TXT = 'The closeness centrality of a node is the inverse of the average shortest path length from the node to all other nodes in the graph.';
        
        ERANGE = 8;
        ERANGE_NAME = 'ERange';
        ERANGE_TXT = 'Calculate the edges which significantly reduce the characteristic path length in the network';
        
        FINDWALK = 9;
        FINDWALK_NAME = 'Find Walk';
        FINDWALK_TXT = 'Walks are sequences of linked nodes, that may visit a single node more than once';
        
        MODULARITY = 10;
        MODULARITY_NAME = 'Modularity';
        MODULARITY_TXT = 'Calculate the edges which significantly reduce the characteristic path length in the network';
        
        RICHCLUB = 11;
        RICHCLUB_NAME = 'Rich Club';
        RICHCLUB_TXT = 'Calculate the edges which significantly reduce the characteristic path length in the network';
        
        Name = { ...
            Project.DEGREE_NAME ...
            Project.PATHLENGTH_NAME ...
            Project.CLUSTERCOEFFICIENT_NAME ...
            Project.RADIUS_NAME ...
            Project.DIAMETER_NAME ...
            Project.BETWEENNESS_NAME ...
            Project.CLOSENESS_NAME ...
            Project.ERANGE_NAME ...
            Project.FINDWALK_NAME ...
            Project.MODULARITY_NAME ...
            Project.RICHCLUB_NAME ...
            
        };

        TXT = { ...
            Project.DEGREE_TXT ...
            Project.PATHLENGTH_TXT ...
            Project.CLUSTERCOEFFICIENT_TXT ...
            Project.RADIUS_TXT ...
            Project.DIAMETER_TXT ...
            Project.BETWEENNESS_TXT ...
            Project.CLOSENESS_TXT ...
            Project.ERANGE_TXT ...
            Project.FINDWALK_TXT ...
            Project.MODULARITY_TXT ...
            Project.RICHCLUB_TXT ...
        };
    end
    
    properties (GetAccess = public, SetAccess = protected)    
        Files       % File of Subject
        L           % Length of Files
        
        Path        % Path of Files
        FileNames   % Name of Files
        Location    % Path to file with the location of the nodes in plot
        Indices     % Indices
        
        Degree      % Average degree of the subjects
        PathLength  % Average path length of the subjects
        ClusterCoefficient % Average cluster coefficient of the subjects
        Radius
        Diameter
        Closeness
        Betweenness
        ERange
        FindWalk
        Modularity
        RichClub

    end
    methods 
        function p = Project(F, Path, FileNames, Location, Indices)

            p.Files = F;
            p.L = length(p.Files);
            
            p.Path = Path;
            p.FileNames = FileNames;
            p.Location = Location;
            p.Indices = Indices;
            
     
        end
    end
    methods
        function runProject(p)
            p.Degree = mean([p.Files.MeanDegree]);
            p.PathLength = mean([p.Files.MeanPathLength]);
            p.ClusterCoefficient = mean([p.Files.MeanClusterCoefficient]);
            p.Radius = mean([p.Files.MeanRadius]);
            p.Diameter = mean([p.Files.MeanDiameter]);
            p.Closeness = mean([p.Files.MeanCloseness]);
            p.Betweenness = mean([p.Files.MeanBetweenness]);
            p.ERange = mean([p.Files.MeanERange]);
            p.FindWalk =  mean([p.Files.MeanFindWalk]);
            p.Modularity = mean([p.Files.MeanModularity]);
            p.RichClub = mean([p.Files.MeanRichClub]);
        end
        
        function write(p)
            
            F = p.Files;
            
            for i=1:length(p.Files)
                relatory(i,:) = {F(i).File ...
                    F(i).MeanDegree F(i).CVDegree  ...
                    F(i).MeanPathLength  F(i).CVPathLength  ...
                    F(i).MeanClusterCoefficient  F(i).CVClusterCoefficient  ...
                    F(i).MeanEdges F(i).CVEdges ...
                    };
            end
             writeFile(relatory,{'Name','Degree', 'CV_Degree', 'PL','CV_PL', ...
                 'CC', 'CV_CC','Edges', 'CV_Edges'},[p.Path '\' 'Relatory'], '.dat');
        end
    end
    

end