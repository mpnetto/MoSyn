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
        
        RADIUS = 4;
        RADIUS_NAME = 'Radius';
        RADIUS_TXT = 'The radius is the minimum eccentricity.';
        RADIUS_DEFAULT = 0;
        
        DIAMETER = 5;
        DIAMETER_NAME = 'Diameter';
        DIAMETER_TXT = 'The diameter is the maximum eccentricity.'
        DIAMETER_DEFAULT = 0;
        
        BETWEENNESS = 6;
        BETWEENNESS_NAME = 'Betweenness';
        BETWEENNESS_TXT = 'Node betweenness centrality of a node is the fraction of all shortest paths in the graph that contain a given node. Nodes with high values of betweenness centrality participate in a large number of shortest paths.';
        BETWEENNESS_DEFAULT = 0;
        
        CLOSENESS = 7;
        CLOSENESS_NAME = 'Closeness';
        CLOSENESS_TXT = 'The closeness centrality of a node is the inverse of the average shortest path length from the node to all other nodes in the graph.';
        CLOSENESS_DEFAULT = 0;
        
        ERANGE = 8;
        ERANGE_NAME = 'ERange';
        ERANGE_TXT = 'Calculate the edges which significantly reduce the characteristic path length in the network';
        ERANGE_DEFAULT = 0;
        
        FINDWALK = 9;
        FINDWALK_NAME = 'Find Walk';
        FINDWALK_TXT = 'Walks are sequences of linked nodes, that may visit a single node more than once';
        FINDWALK_DEFAULT = 0;
        
        MODULARITY = 10;
        MODULARITY_NAME = 'Modularity';
        MODULARITY_TXT = 'Calculate the edges which significantly reduce the characteristic path length in the network';
        MODULARITY_DEFAULT = 0;
        
        RICHCLUB = 11;
        RICHCLUB_NAME = 'Rich Club';
        RICHCLUB_TXT = 'Calculate the edges which significantly reduce the characteristic path length in the network';
        RICHCLUB_DEFAULT = 0;
        
        Name = { ...
            TVG.DEGREE_NAME ...
            TVG.PATHLENGTH_NAME ...
            TVG.CLUSTERCOEFFICIENT_NAME ...
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
        T           % tensor connection matrix
        T_W         % tensor connection weighted matrix
        T_D_W       % tensor connection oriented weighted matrix
        L           % length of TVG

        File        % name of the file
        Path        % path for the fila
        Labels      % labels of the nodes
        Indices     % Indeces to calculate
        Location
        
        InitialTime % initial time of the time series
        FinalTime   % final time of the time series
        NumNodes    % amount of nodes
        SlidWindow  % windows size
        TaoMin      % minimal size of the delay time
        TaoMax      % maximal size of the delay time
        Threshold   % cut limit of the synchronization
        OutTvg      % Flag to calc or not all TVG edges
        TSCLim      % Limt fot TSC calc
        
        Degree      %
        Weight      %
        Hub         %
        Edges       %
        tsc
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
        
        MicroStates
        
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
        function t = TVG(File, Path, T, T_W,T_D_W, Labels, Indices, InitialTime, FinalTime, NumNodes, SlidWindow, TaoMin, TaoMax, Threshold,tscLim,OutTvg,Location)
            if nargin > 0
                t.T = T;
                t.T_W = T_W;
                t.T_D_W = T_D_W;
                t.L = length(t.T(1,1,:));
                t.File = File;
                t.Path = Path;
                t.Indices = Indices;
                t.InitialTime = InitialTime;
                t.FinalTime = FinalTime;
                t.NumNodes = NumNodes;
                t.SlidWindow = SlidWindow;
                t.TaoMin = TaoMin;
                t.TaoMax = TaoMax;
                t.Threshold = Threshold;
                t.Labels = Labels;
                t.OutTvg = OutTvg;
                t.TSCLim=tscLim;
                t.Location = Location;
            end
            
        end
        function error = runTVG(tvg)
            
            error = 0;
            T = tvg.T;
            T_W = tvg.T_W;
            
            locationNodes = readLocation(tvg);
            
            if tvg.OutTvg
                T_D_W = tvg.T_D_W;
                Tau = cell2table(cell(0,5), 'VariableNames', {'time','Source','Target','Tau','Distance'});
            end
            
            indices = tvg.Indices;
            HomoEle=zeros(tvg.NumNodes,2);

            tStart = tic;
            hbar = parfor_progressbar(tvg.L,'Computing Network');  %create the progress bar
            %multiWaitbar('Network',0, 'Color', 'b', 'CanCancel','on', 'CancelFcn', @(a,b) disp( ['Cancel ',a] ));
            if tvg.OutTvg
                parfor j=1:tvg.L  %parfor
                    hbar.iterate(1);
                    
                    M = T(:,:,j);
                    MicroStates(j,:) = M(:);
                    
                    A = T_D_W(:,:,j);
                    t = tvg.TAUcalc(double(A), locationNodes,j);
                    if(~isempty(t))
                        Tau = vertcat(Tau, t);
                    end
                    
                    g = GraphM(M);
                    
                    if find(strcmp(indices, tvg.DEGREE_NAME))
                        [Degree(:,j),MeanDegree(:,j), ...
                            StdDegree(:,j)]  = g.binDegree();
                    end
                    
                    if find(strcmp(indices, tvg.CLUSTERCOEFFICIENT_NAME))
                        [Cluster(:,j)] = g.cluster();
                    end
                    if find(strcmp(indices, tvg.PATHLENGTH_NAME))
                        PathLength(:,j)  = g.pathLength();
                    end
                    
                    if find(strcmp(indices, tvg.RADIUS_NAME))
                        Radius(:,j)  = g.radius();
                    end
                    
                    if find(strcmp(indices, tvg.DIAMETER_NAME))
                        Diameter(:,j)  = g.diameter();
                    end
                    
                    if find(strcmp(indices, tvg.CLOSENESS_NAME))
                        Closeness(:,j)  = g.closeness();
                    end
                    
                    if find(strcmp(indices, tvg.BETWEENNESS_NAME))
                        Betweenness(:,j)  = g.betweenness();
                    end
                    
                    if find(strcmp(indices, tvg.ERANGE_NAME))
                        [~,~,~,ERange(:,j)]  = g.erange();
                    end
                    
                    if find(strcmp(indices, tvg.FINDWALK_NAME))
                        [~,FindWalk(:,j)]  = g.findwalks();
                    end
                    
                    if find(strcmp(indices, tvg.MODULARITY_NAME))
                        [~,Modularity(:,j)]  = g.modularity_und(1);
                    end
                    
                    if find(strcmp(indices, tvg.RICHCLUB_NAME))
                        [~,aux]  = g.rich_club_bd();
                        RichClub(:,j) = mean(aux);
                    end
                    
                    h = 0; %g.eIndex(tvg.Labels);
                    HomoEle = HomoEle + h;
                    
                    HubDeg(:,j)= g.hub()';
                    Edges(:,j) = g.edges();
                    
                end
            else
                parfor j=1:tvg.L  %parfor
                    hbar.iterate(1);
                    M = T(:,:,j);
                    MicroStates(j,:) = M(:);
                    g = GraphM(M);
                    
                    if find(strcmp(indices, tvg.DEGREE_NAME))
                        [Degree(:,j),MeanDegree(:,j), ...
                            StdDegree(:,j)]  = g.binDegree();
                    end
                    
                    if find(strcmp(indices, tvg.CLUSTERCOEFFICIENT_NAME))
                        [Cluster(:,j)] = g.cluster();
                    end
                    if find(strcmp(indices, tvg.PATHLENGTH_NAME))
                        PathLength(:,j)  = g.pathLength();
                    end
                    
                    if find(strcmp(indices, tvg.RADIUS_NAME))
                        Radius(:,j)  = g.radius();
                    end
                    
                    if find(strcmp(indices, tvg.DIAMETER_NAME))
                        Diameter(:,j)  = g.diameter();
                    end
                    
                    if find(strcmp(indices, tvg.CLOSENESS_NAME))
                        Closeness(:,j)  = g.closeness();
                    end
                    
                    if find(strcmp(indices, tvg.BETWEENNESS_NAME))
                        Betweenness(:,j)  = g.betweenness();
                    end
                    
                    if find(strcmp(indices, tvg.ERANGE_NAME))
                        [~,~,~,ERange(:,j)]  = g.erange();
                    end
                    
                    if find(strcmp(indices, tvg.FINDWALK_NAME))
                        [~,FindWalk(:,j)]  = g.findwalks();
                    end
                    
                    if find(strcmp(indices, tvg.MODULARITY_NAME))
                        [~,Modularity(:,j)]  = g.modularity_und(1);
                    end
                    
                    if find(strcmp(indices, tvg.RICHCLUB_NAME))
                        [~,aux]  = g.rich_club_bd();
                        RichClub(:,j) = mean(aux);
                    end
                    
                    h = 0; % g.eIndex(tvg.Labels);
                    HomoEle = HomoEle + h;
                    
                    HubDeg(:,j)= g.hub()';
                    Edges(:,j) = g.edges();
                    
                end
                
            end

            close(hbar);

            tvg.Hub = HubDeg;
            tvg.Edges = Edges;
            if (tvg.TSCLim>0)
                tvg.tsc = tsc(T, tvg.TSCLim);
            else
                tvg.tsc = -1;   %% not calculated
            end
            
            if find(strcmp(indices, tvg.DEGREE_NAME))
                tvg.Degree = Degree;
                tvg.MeanFrameDegree = mean(Degree)';
                tvg.MeanDegree = mean(tvg.MeanFrameDegree);
                tvg.StdDegree = std(tvg.MeanFrameDegree);
                tvg.CVDegree =  tvg.StdDegree / tvg.MeanDegree;
            end
            
            if find(strcmp(indices, tvg.PATHLENGTH_NAME))
                tvg.PathLength = PathLength;
                tvg.MeanFramePathLength = mean(PathLength)';
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
            
            if find(strcmp(indices, tvg.RADIUS_NAME))
                tvg.Radius = Radius;
                tvg.MeanFrameRadius = tvg.Radius;
                tvg.MeanRadius = mean(tvg.Radius);
                tvg.StdRadius = std(tvg.Radius);
                tvg.CVRadius = tvg.StdRadius / tvg.MeanRadius;
            end

            if find(strcmp(indices, tvg.DIAMETER_NAME))
                tvg.Diameter = Diameter;
                tvg.MeanFrameDiameter = tvg.Diameter;
                tvg.MeanDiameter = mean(tvg.Diameter);
                tvg.StdDiameter = std(tvg.Diameter);
                tvg.CVDiameter = tvg.StdDiameter / tvg.MeanDiameter;
            end

            if find(strcmp(indices, tvg.CLOSENESS_NAME))
                tvg.Closeness = Closeness;
                tvg.MeanFrameCloseness = mean(tvg.Closeness)';
                tvg.MeanCloseness = mean(tvg.MeanFrameCloseness);
                tvg.StdCloseness = std(tvg.MeanFrameCloseness);
                tvg.CVCloseness = tvg.StdCloseness / tvg.MeanCloseness;
            end

            if find(strcmp(indices, tvg.BETWEENNESS_NAME))
                tvg.Betweenness = Betweenness;
                tvg.MeanFrameBetweenness = mean(tvg.Betweenness)';
                tvg.MeanBetweenness = mean(tvg.MeanFrameBetweenness);
                tvg.StdBetweenness = std(tvg.MeanFrameBetweenness);
                tvg.CVBetweenness = tvg.StdBetweenness / tvg.MeanBetweenness;
            end
            
            if find(strcmp(indices, tvg.ERANGE_NAME))
                tvg.ERange = ERange;
                tvg.MeanFrameERange = tvg.ERange;
                tvg.MeanERange = mean(tvg.MeanFrameERange);
                tvg.StdERange = std(tvg.MeanFrameERange);
                tvg.CVERange = tvg.StdERange / tvg.MeanERange;
            end
            
            if find(strcmp(indices, tvg.FINDWALK_NAME))
                tvg.FindWalk = FindWalk;
                tvg.MeanFrameFindWalk = tvg.FindWalk;
                tvg.MeanFindWalk = mean(tvg.MeanFrameFindWalk);
                tvg.StdFindWalk = std(tvg.MeanFrameFindWalk);
                tvg.CVFindWalk = tvg.StdFindWalk / tvg.MeanFindWalk;
            end

            if find(strcmp(indices, tvg.MODULARITY_NAME))
                tvg.Modularity = Modularity;
                tvg.MeanFrameModularity =  tvg.Modularity;
                tvg.MeanModularity = mean(tvg.MeanFrameModularity);
                tvg.StdModularity = std(tvg.MeanFrameModularity);
                tvg.CVModularity = tvg.StdModularity / tvg.MeanModularity;
            end

            if find(strcmp(indices, tvg.RICHCLUB_NAME))
                tvg.RichClub = RichClub;
                tvg.MeanFrameRichClub =tvg.RichClub;
                tvg.MeanRichClub = mean(tvg.MeanFrameRichClub);
                tvg.StdRichClub = std(tvg.MeanFrameRichClub);
                tvg.CVRichClub = tvg.StdRichClub / tvg.MeanRichClub;
            end
            
            %A = logical(MicroStates);
            %B = bi2de(A);
            %tvg.MicroStates = B;
            
            tvg.MeanEdges = mean(tvg.Edges);
            tvg.StdEdges = std(tvg.Edges);
            tvg.CVEdges = tvg.StdEdges / tvg.MeanEdges;
            
            tvg.REA = REA(sum(T,3), sum(T_W,3), tvg.File, tvg.Path, tvg.Labels);
            tvg.REA.runREA();
            tvg.REA.write();
            
            Degree_In = squeeze(sum(T_W,1));
            HubDeg_In = Degree_In > (mean(Degree_In)+ 2*std(Degree_In,1));
            Degree_Out = squeeze(sum(T_W,2));
            HubDeg_Out = Degree_Out > (mean(Degree_Out)+ 2*std(Degree_Out,1));
            
            tvg.HubSim = sum(HubDeg,2);
            tvg.Hubs = tvg.HubSim./sum(tvg.HubSim,1);
            tvg.HubIn = sum(HubDeg_In,2);
            tvg.HubsIn = tvg.HubIn./sum(tvg.HubIn,1);
            tvg.HubOut = sum(HubDeg_Out,2);
            tvg.HubsOut = tvg.HubOut./sum(tvg.HubOut,1);
            if tvg.OutTvg
                Freq=histcounts(Tau.Tau,tvg.TaoMax-tvg.TaoMin+1)';
                Freq(1)=Freq(1)/2; % disregard double count in simetric edges
                ctau=(tvg.TaoMin+1:tvg.TaoMax+1)';
                Ht=[array2table(ctau) array2table(Freq)];
                NameSync = 'MoS';
                writetable(Ht, [tvg.Path '\' tvg.File '_TAU_Hist_' NameSync  '.txt'], 'Delimiter','\t');
                writetable(Tau, [tvg.Path '\' tvg.File '_TAU_' NameSync  '.txt'], 'Delimiter','\t');
            end

        end
        
        function write(t)
            NameSync = 'MoS';
            L = t.L - 1;
            indices = t.Indices;
            
            DistHubs = horzcat(t.Labels,num2cell(t.HubSim), num2cell(t.HubIn), num2cell(t.HubOut), num2cell(t.Hubs), num2cell(t.HubsIn), num2cell(t.HubsOut));  
            TimeSerie=horzcat((0:1:L)',t.Edges',t.MeanFrameClusterCoefficient);
            TimeStat=horzcat(t.MeanEdges,t.StdEdges,t.MeanClusterCoefficient,t.StdClusterCoefficient);

            if find(strcmp(indices, t.DEGREE_NAME))
                writeFile(horzcat((0:1:L)', t.MeanFrameDegree),{'Frame','MeanDegree'},[t.Path '\' t.File '_Degree'], '.txt');
            end
            if find(strcmp(indices, t.PATHLENGTH_NAME))
                writeFile(horzcat((0:1:L)', t.MeanFramePathLength),{'Frame','PL'},[t.Path '\' t.File '_PathLength'], '.txt');
            end
            if find(strcmp(indices, t.CLUSTERCOEFFICIENT_NAME))
                writeFile(horzcat((0:1:L)', t.MeanFrameClusterCoefficient),{'Frame','CC'},[t.Path '\' t.File '_CC'], '.txt'); 
                writeFile( ...
                    TimeSerie, ...
                    {'Time' 'Edges' 'CC'}, ...
                    [t.Path '\' t.File '_iTime_' NameSync], ...
                    '.txt');
                 writeFile( ...
                    TimeStat, ...
                    {'edgMed' 'edgStd' 'cluMed', 'cluStd'}, ...
                    [t.Path '\' t.File '_sTime_' NameSync], ...
                    '.txt');
            else
                writeFile( ...
                    TimeSerie, ...
                    {'Time' 'Edges'}, ...
                    [t.Path '\' t.File '_iTime_' NameSync], ...
                    '.txt');
                writeFile( ...
                    TimeStat, ...
                    {'edgMed' 'edgStd'}, ...
                    [t.Path '\' t.File '_sTime_' NameSync], ...
                    '.txt');
            
            end
            
            
            if(t.tsc>=0)
                writeFile(t.tsc,{'TSC'},[t.Path '\' t.File '_TSC'], '.txt');
            end
            
            
             writeFile( ...
                DistHubs, ...
                {'Label', 'HubSim', 'HubIn', 'HubOut',	'Hubs',	'HubsIn', 'HubsOut'}, ...
                [t.Path '\' t.File '_Hubs_' NameSync], ...
                '.txt');
          
           
        end
        function tau = TAUcalc(tvg, A, locationNodes, Time)
            
            labels = locationNodes.Properties.RowNames;
            
            G = digraph(A,labels);
            e = splitvars(G.Edges);
            
            for z = 1:height(e)
                row = e(z,:);
                node1 = table2cell(row(:,1));
                node2 = table2cell(row(:,2));
                
                coordNode1 = table2array(locationNodes(node1,:));
                coordNode2 = table2array(locationNodes(node2,:));
                
                Distance(z,:) = norm(coordNode2 - coordNode1);
            end
            if (height(e)>0)
                ti(1:height(e))=Time;
                tau = [array2table(ti') e array2table(Distance)];
                tau.Properties.VariableNames = {'time','Source','Target','Tau','Distance'};
            else
                tau=[];
            end
            
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