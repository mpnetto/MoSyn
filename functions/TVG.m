classdef TVG < handle
   
    properties (Constant)
            
        DEGREE = 1;
        DEGREE_NAME = 'Degree';
        DEGREE_TXT = 'The degree of a node is the number of edges connceted to the node.';
        DEGREE_DEFAULT = 1;
        
        PATHLENGTH = 2;
        PATHLENGTH_NAME = 'Path Length';
        PATHLENGTH_TXT = 'The path length of a node is the average minimum path between on node to every other node.';
        PATHLENGTH_DEFAULT = 1;
        
        CLUSTERCOEFFICIENT = 3;
        CLUSTERCOEFFICIENT_NAME = 'Cluster Coefficient';
        CLUSTERCOEFFICIENT_TXT = 'The clustering coefficient of a graph is the average of the clustering coefficients of its nodes.';
        CLUSTERCOEFFICIENT_DEFAULT = 1;

        TAU = 4;
        TAU_NAME = 'Tau';
        TAU_TXT = 'Temporary tau explanation';
        TAU_DEFAULT = 0;
        
        RADIUS = 5;
        RADIUS_NAME = 'Radius';
        RADIUS_TXT = 'The radius is the minimum eccentricity.';
        RADIUS_DEFAULT = 0;
        
        DIAMETER = 6;
        DIAMETER_NAME = 'Diameter';
        DIAMETER_TXT = 'The diameter is the maximum eccentricity.'
        DIAMETER_DEFAULT = 0;
        
        BETWEENNESS = 7;
        BETWEENNESS_NAME = 'Betweenness';
        BETWEENNESS_TXT = 'Node betweenness centrality of a node is the fraction of all shortest paths in the graph that contain a given node. Nodes with high values of betweenness centrality participate in a large number of shortest paths.';
        BETWEENNESS_DEFAULT = 0;
        
        CLOSENESS = 8;
        CLOSENESS_NAME = 'Closeness';
        CLOSENESS_TXT = 'The closeness centrality of a node is the inverse of the average shortest path length from the node to all other nodes in the graph.';
        CLOSENESS_DEFAULT = 0;
        
        ERANGE = 9;
        ERANGE_NAME = 'ERange';
        ERANGE_TXT = 'Calculate the edges which significantly reduce the characteristic path length in the network';
        ERANGE_DEFAULT = 0;
        
        FINDWALK = 10;
        FINDWALK_NAME = 'Find Walk';
        FINDWALK_TXT = 'Walks are sequences of linked nodes, that may visit a single node more than once';
        FINDWALK_DEFAULT = 0;
        
        MODULARITY = 11;
        MODULARITY_NAME = 'Modularity';
        MODULARITY_TXT = 'Calculate the edges which significantly reduce the characteristic path length in the network';
        MODULARITY_DEFAULT = 0;
        
        RICHCLUB = 12;
        RICHCLUB_NAME = 'Rich Club';
        RICHCLUB_TXT = 'Calculate the edges which significantly reduce the characteristic path length in the network';
        RICHCLUB_DEFAULT = 0;
        
        Name = { ...
            TVG.DEGREE_NAME ...
            TVG.PATHLENGTH_NAME ...
            TVG.CLUSTERCOEFFICIENT_NAME ...
            TVG.TAU_NAME ...
            TVG.RADIUS_NAME ...
            TVG.DIAMETER_NAME ...
            TVG.BETWEENNESS_NAME ...
            TVG.CLOSENESS_NAME ...
            TVG.ERANGE_NAME ...
            TVG.FINDWALK_NAME ...
            TVG.MODULARITY_NAME ...
            TVG.RICHCLUB_NAME ...
            
        };

        TXT = { ...
            TVG.DEGREE_TXT ...
            TVG.PATHLENGTH_TXT ...
            TVG.CLUSTERCOEFFICIENT_TXT ...
            TVG.TAU_NAME ...
            TVG.RADIUS_TXT ...
            TVG.DIAMETER_TXT ...
            TVG.BETWEENNESS_TXT ...
            TVG.CLOSENESS_TXT ...
            TVG.ERANGE_TXT ...
            TVG.FINDWALK_TXT ...
            TVG.MODULARITY_TXT ...
            TVG.RICHCLUB_TXT ...
        };
    
        CHECK = { ...
            TVG.DEGREE_DEFAULT ...
            TVG.PATHLENGTH_DEFAULT ...
            TVG.CLUSTERCOEFFICIENT_DEFAULT ...
            TVG.TAU_DEFAULT ...
            TVG.RADIUS_DEFAULT ...
            TVG.DIAMETER_DEFAULT ...
            TVG.BETWEENNESS_DEFAULT ...
            TVG.CLOSENESS_DEFAULT ...
            TVG.ERANGE_DEFAULT ...
            TVG.FINDWALK_DEFAULT ...
            TVG.MODULARITY_DEFAULT ...
            TVG.RICHCLUB_DEFAULT ...
        };
    end
    
    properties (GetAccess = public, SetAccess = protected)
        TVGArray
        WeightedTVG
        OrientedWeightedTVG
        Length

        Filename
        FilePath
        NodesLabels
        Indices
        NodesLocation
        
        InitialTime
        FinalTime
        NumNodes
        SlidWindow
        TaoMin
        TaoMax
        Threshold
        TSCLimit
        
        Degree
        Weight
        Hub
        Edges
        TSC
        MeanEdges   
        StdEdges
        CVEdges
        PathLength
        ClusterCoefficient
        HubSim
        Hubs
        HubIn
        HubsIn
        HubOut
        HubsOut
        
        MeanFrameDegree
        MeanFramePathLength
        MeanFrameClusterCoefficient
        
        MeanDegree
        StdDegree
        CVDegree
        
        MeanPathLength
        StdPathLength
        CVPathLength
        
        MeanClusterCoefficient
        StdClusterCoefficient
        CVClusterCoefficient

        Tau
        TauDist
        MeanTau
        MeanFrameTau
        StdTau
        CVTau
        
        Radius
        MeanRadius
        MeanFrameRadius
        StdRadius
        CVRadius

        Diameter
        MeanDiameter
        MeanFrameDiameter
        StdDiameter
        CVDiameter

        Closeness
        MeanFrameCloseness
        MeanCloseness
        StdCloseness
        CVCloseness

        Betweenness
        MeanFrameBetweenness
        MeanBetweenness
        StdBetweenness
        CVBetweenness
        
        ERange
        MeanFrameERange
        MeanERange
        StdERange
        CVERange
        
        FindWalk
        MeanFrameFindWalk
        MeanFindWalk
        StdFindWalk
        CVFindWalk
        
        Modularity
        MeanFrameModularity
        MeanModularity
        StdModularity
        CVModularity
        
        RichClub
        MeanFrameRichClub
        MeanRichClub
        StdRichClub
        CVRichClub
  
        REA
 
    end
    methods

        % Get and save the TVG parameters
        function t = TVG(File, Path, T, T_W,T_D_W, Labels, Indices, InitialTime, FinalTime, NumNodes, SlidWindow, TaoMin, TaoMax, Threshold, tscLim, Location)
            if nargin > 0
                t.TVGArray = T;
                t.WeightedTVG = T_W;
                t.OrientedWeightedTVG = T_D_W;
                t.Length = length(t.TVGArray(1,1,:));
                t.Filename = File;
                t.FilePath = Path;
                t.Indices = Indices;
                t.InitialTime = InitialTime;
                t.FinalTime = FinalTime;
                t.NumNodes = NumNodes;
                t.SlidWindow = SlidWindow;
                t.TaoMin = TaoMin;
                t.TaoMax = TaoMax;
                t.Threshold = Threshold;
                t.NodesLabels = Labels;
                t.TSCLimit = tscLim;
                t.NodesLocation = Location;
            end
            
        end
        
        function runTVG(tvg)
            
            tvgArray = tvg.TVGArray;
            weightedTvg = tvg.WeightedTVG;
            orientedWeightedTVG = tvg.OrientedWeightedTVG;
            indices = tvg.Indices;
            locationNodes = readLocation(tvg);
            
            tau = [] ; % cell2table(cell(0,5), 'VariableNames', {'time','Source','Target','Tau','Distance'});
            tauDistance = [] ;
            for j=1:tvg.Length
                
                M = tvgArray(:,:,j);

                g = GraphM(M);
                
                if find(strcmp(indices, tvg.DEGREE_NAME))
                    [degree(:,j), ~, ~] = g.binDegree();
                end
                
                if find(strcmp(indices, tvg.CLUSTERCOEFFICIENT_NAME))
                    [Cluster(:,j)] = g.cluster();
                end

                if find(strcmp(indices, tvg.PATHLENGTH_NAME))
                    pathLength(:,j) = g.pathLength();
                end

                if find(strcmp(indices, tvg.TAU_NAME))
                    a = orientedWeightedTVG(:,:,j);
                    ag = GraphM(a, locationNodes, j);

                    [tauValues, tauDist] = ag.calcTau();
                    tau = vertcat(tau, tauValues);
                    tauDistance = vertcat(tauDistance, tauDist);

                end
                
                if find(strcmp(indices, tvg.RADIUS_NAME))
                    radius(:,j) = g.radius();
                end
                
                if find(strcmp(indices, tvg.DIAMETER_NAME))
                    diameter(:,j) = g.diameter();
                end
                
                if find(strcmp(indices, tvg.CLOSENESS_NAME))
                    closeness(:,j) = g.closeness();
                end
                
                if find(strcmp(indices, tvg.BETWEENNESS_NAME))
                    betweenness(:,j) = g.betweenness();
                end
                
                if find(strcmp(indices, tvg.ERANGE_NAME))
                    [~,~,~,eRange(:,j)] = g.erange();
                end
                
                if find(strcmp(indices, tvg.FINDWALK_NAME))
                    [~,findWalk(:,j)] = g.findwalks();
                end
                
                if find(strcmp(indices, tvg.MODULARITY_NAME))
                    [~,modularity(:,j)] = g.modularity_und(1);
                end
                
                if find(strcmp(indices, tvg.RICHCLUB_NAME))
                    [~,aux]  = g.rich_club_bd();
                    richClub(:,j) = mean(aux);
                end
                
                hubDeg(:,j)= g.hub()';
                edges(:,j) = g.edges();
                
            end

            tvg.Hub = hubDeg;
            tvg.Edges = edges;
            
            if (tvg.TSCLimit > 0)
                tvg.TSC = tsc(tvgArray, tvg.TSCLimit);
            end
            
            if find(strcmp(indices, tvg.DEGREE_NAME))
                tvg.Degree = degree;
                tvg.MeanFrameDegree = mean(degree)';
                tvg.MeanDegree = mean(tvg.MeanFrameDegree);
                tvg.StdDegree = std(tvg.MeanFrameDegree);
                tvg.CVDegree =  tvg.StdDegree / tvg.MeanDegree;
            end
            
            if find(strcmp(indices, tvg.PATHLENGTH_NAME))
                tvg.PathLength = pathLength;
                tvg.MeanFramePathLength = mean(pathLength)';
                tvg.MeanPathLength = mean(tvg.MeanFramePathLength);
                tvg.StdPathLength = std(tvg.MeanFramePathLength);
                tvg.CVPathLength = tvg.StdPathLength /tvg.MeanPathLength;
            end
            
            if find(strcmp(indices, tvg.CLUSTERCOEFFICIENT_NAME))
                tvg.ClusterCoefficient = Cluster;
                tvg.MeanFrameClusterCoefficient = mean(Cluster)';
                tvg.MeanClusterCoefficient = mean(tvg.MeanFrameClusterCoefficient);
                tvg.StdClusterCoefficient = std(tvg.MeanFrameClusterCoefficient);
                tvg.CVClusterCoefficient = tvg.StdClusterCoefficient / tvg.MeanClusterCoefficient;
            end

             if find(strcmp(indices, tvg.TAU_NAME))
                tvg.Tau = tau;
                tvg.TauDist = tauDistance;
%                 tvg.MeanFrameTau = mean(tau)';
%                 tvg.MeanTau = mean(tvg.MeanFrameTau);
%                 tvg.StdTau = std(tvg.MeanTau);
%                 tvg.CVTau = tvg.StdTau / tvg.MeanTau;
            end
            
            if find(strcmp(indices, tvg.RADIUS_NAME))
                tvg.Radius = radius;
                tvg.MeanFrameRadius = tvg.Radius;
                tvg.MeanRadius = mean(tvg.Radius);
                tvg.StdRadius = std(tvg.Radius);
                tvg.CVRadius = tvg.StdRadius / tvg.MeanRadius;
            end

            if find(strcmp(indices, tvg.DIAMETER_NAME))
                tvg.Diameter = diameter;
                tvg.MeanFrameDiameter = tvg.Diameter;
                tvg.MeanDiameter = mean(tvg.Diameter);
                tvg.StdDiameter = std(tvg.Diameter);
                tvg.CVDiameter = tvg.StdDiameter / tvg.MeanDiameter;
            end

            if find(strcmp(indices, tvg.CLOSENESS_NAME))
                tvg.Closeness = closeness;
                tvg.MeanFrameCloseness = mean(tvg.Closeness)';
                tvg.MeanCloseness = mean(tvg.MeanFrameCloseness);
                tvg.StdCloseness = std(tvg.MeanFrameCloseness);
                tvg.CVCloseness = tvg.StdCloseness / tvg.MeanCloseness;
            end

            if find(strcmp(indices, tvg.BETWEENNESS_NAME))
                tvg.Betweenness = betweenness;
                tvg.MeanFrameBetweenness = mean(tvg.Betweenness)';
                tvg.MeanBetweenness = mean(tvg.MeanFrameBetweenness);
                tvg.StdBetweenness = std(tvg.MeanFrameBetweenness);
                tvg.CVBetweenness = tvg.StdBetweenness / tvg.MeanBetweenness;
            end
            
            if find(strcmp(indices, tvg.ERANGE_NAME))
                tvg.ERange = eRange;
                tvg.MeanFrameERange = tvg.ERange;
                tvg.MeanERange = mean(tvg.MeanFrameERange);
                tvg.StdERange = std(tvg.MeanFrameERange);
                tvg.CVERange = tvg.StdERange / tvg.MeanERange;
            end
            
            if find(strcmp(indices, tvg.FINDWALK_NAME))
                tvg.FindWalk = findWalk;
                tvg.MeanFrameFindWalk = tvg.FindWalk;
                tvg.MeanFindWalk = mean(tvg.MeanFrameFindWalk);
                tvg.StdFindWalk = std(tvg.MeanFrameFindWalk);
                tvg.CVFindWalk = tvg.StdFindWalk / tvg.MeanFindWalk;
            end

            if find(strcmp(indices, tvg.MODULARITY_NAME))
                tvg.Modularity = modularity;
                tvg.MeanFrameModularity =  tvg.Modularity;
                tvg.MeanModularity = mean(tvg.MeanFrameModularity);
                tvg.StdModularity = std(tvg.MeanFrameModularity);
                tvg.CVModularity = tvg.StdModularity / tvg.MeanModularity;
            end

            if find(strcmp(indices, tvg.RICHCLUB_NAME))
                tvg.RichClub = richClub;
                tvg.MeanFrameRichClub =tvg.RichClub;
                tvg.MeanRichClub = mean(tvg.MeanFrameRichClub);
                tvg.StdRichClub = std(tvg.MeanFrameRichClub);
                tvg.CVRichClub = tvg.StdRichClub / tvg.MeanRichClub;
            end
           
            tvg.MeanEdges = mean(tvg.Edges);
            tvg.StdEdges = std(tvg.Edges);
            tvg.CVEdges = tvg.StdEdges / tvg.MeanEdges;
            
            tvg.REA = REA(sum(tvgArray,3), sum(weightedTvg,3), tvg.Filename, tvg.FilePath, tvg.NodesLabels);
            tvg.REA.runREA();
            tvg.REA.write();
            
            Degree_In = squeeze(sum(weightedTvg,1));
            HubDeg_In = Degree_In > (mean(Degree_In)+ 2*std(Degree_In,1));
            Degree_Out = squeeze(sum(weightedTvg,2));
            HubDeg_Out = Degree_Out > (mean(Degree_Out)+ 2*std(Degree_Out,1));
            
            tvg.HubSim = sum(hubDeg,2);
            tvg.Hubs = tvg.HubSim./sum(tvg.HubSim,1);
            tvg.HubIn = sum(HubDeg_In,2);
            tvg.HubsIn = tvg.HubIn./sum(tvg.HubIn,1);
            tvg.HubOut = sum(HubDeg_Out,2);
            tvg.HubsOut = tvg.HubOut./sum(tvg.HubOut,1);

        end
        
        function write(tvg)
            NameSync = 'MoS';
            length = tvg.Length - 1;
            indices = tvg.Indices;
            
            DistHubs = horzcat(tvg.NodesLabels, num2cell(tvg.HubSim), num2cell(tvg.HubIn), num2cell(tvg.HubOut), num2cell(tvg.Hubs), num2cell(tvg.HubsIn), num2cell(tvg.HubsOut));  
            TimeSerie = horzcat((0:1:length)',tvg.Edges',tvg.MeanFrameClusterCoefficient);
            TimeStat = horzcat(tvg.MeanEdges,tvg.StdEdges,tvg.MeanClusterCoefficient,tvg.StdClusterCoefficient);

            if find(strcmp(indices, tvg.DEGREE_NAME))
                writeFile(horzcat((0:1:length)', tvg.MeanFrameDegree),{'Frame','MeanDegree'},[tvg.FilePath '\' tvg.Filename '_Degree'], '.txt');
            end

            if find(strcmp(indices, tvg.PATHLENGTH_NAME))
                writeFile(horzcat((0:1:length)', tvg.MeanFramePathLength),{'Frame','PL'},[tvg.FilePath '\' tvg.Filename '_PathLength'], '.txt');
            end

            if find(strcmp(indices, tvg.CLUSTERCOEFFICIENT_NAME))
                writeFile(horzcat((0:1:length)', tvg.MeanFrameClusterCoefficient),{'Frame','CC'},[tvg.FilePath '\' tvg.Filename '_CC'], '.txt'); 
                writeFile( ...
                    TimeSerie, ...
                    {'Time' 'Edges' 'CC'}, ...
                    [tvg.FilePath '\' tvg.Filename '_iTime_' NameSync], ...
                    '.txt');
                 writeFile( ...
                    TimeStat, ...
                    {'edgMed' 'edgStd' 'cluMed', 'cluStd'}, ...
                    [tvg.FilePath '\' tvg.Filename '_sTime_' NameSync], ...
                    '.txt');
            else
                writeFile( ...
                    TimeSerie, ...
                    {'Time' 'Edges'}, ...
                    [tvg.FilePath '\' tvg.Filename '_iTime_' NameSync], ...
                    '.txt');
                writeFile( ...
                    TimeStat, ...
                    {'edgMed' 'edgStd'}, ...
                    [tvg.FilePath '\' tvg.Filename '_sTime_' NameSync], ...
                    '.txt');
            
            end

            if find(strcmp(indices, tvg.TAU_NAME))
                % cell2table(cell(0,5), 'VariableNames', {'time','Source','Target','Tau','Distance'});
                Freq = histcounts(tvg.Tau(:,4), tvg.TaoMax - tvg.TaoMin+1)';
                Freq(1) = Freq(1) / 2; % disregard double count in simetric edges
                ctau=(tvg.TaoMin+1 : tvg.TaoMax+1)';
                Ht=[array2table(ctau) array2table(Freq)];



                tauTable = [array2table(tvg.Tau, 'VariableNames', {'time','Source','Target','Tau'}) array2table(tvg.TauDist, 'VariableNames', {'Distance'})];

                tauTable.Source = tvg.NodesLabels( tauTable.Source);
                tauTable.Target = tvg.NodesLabels( tauTable.Target);


                writetable(Ht, [tvg.FilePath '\' tvg.Filename '_TAU_Hist_' NameSync  '.txt'], 'Delimiter','\t');
                writetable(tauTable, [tvg.FilePath '\' tvg.Filename '_TAU_' NameSync  '.txt'], 'Delimiter','\t');
            end
            
            
            if(tvg.TSC>=0)
                writeFile(tvg.TSC,{'TSC'},[tvg.FilePath '\' tvg.Filename '_TSC'], '.txt');
            end

             writeFile( ...
                DistHubs, ...
                {'Label', 'HubSim', 'HubIn', 'HubOut',	'Hubs',	'HubsIn', 'HubsOut'}, ...
                [tvg.FilePath '\' tvg.Filename '_Hubs_' NameSync], ...
                '.txt');
          
           
        end

        function clearTVG(tvg)
            tvg.WeightedTVG = [];
            tvg.OrientedWeightedTVG = [];
        end

    end
    
    methods (Static)
        function z = zeros(varargin)
            if(nargin == 0)
                z = TVG;
            elseif any([varargin{:}] <= 0)
                z = TVG.empty(varargin{:});
            else
                z = repmat(TVG,varargin{:});
            end
        end
    end
end