classdef REA < handle
   
    properties (Constant)
            
        DEGREE = 1;
        DEGREE_NAME = 'Degree';
        DEGREE_NODAL = 'false';
        DEGREE_TXT = 'The degree of a REA is the average of the degree of each node in REA.';
        
        WEIGHT = 2;
        WEIGHT_NAME = 'Weight';
        WEIGHT_NODAL = 'false';
        WEIGHT_TXT = 'A construir';

        
        Name = { ...
            REA.DEGREE_NAME ...
            REA.WEIGHT_NAME ...
        };

        NODAL = { ...
            REA.DEGREE_NODAL ...
            REA.WEIGHT_NODAL ...
        };
        TXT = { ...
            REA.DEGREE_TXT ...
            REA.WEIGHT_TXT ...
        };
    end
    
    properties (GetAccess = public, SetAccess = protected)
        A % connection matrix
        W
        N
        
        File
        Path
        Labels
        
        Weight
        
        Degree
        MeanDegree
        StdDegree
        
        Cluster
        
        WeiDeg 
        WeiDegIn
        WeiDegOut
        
        InsideCon
        OutsideCon
        Eindex
        
        Lindex
    end
    methods 
        function r = REA(A, W, File, Path, Labels)
            r.A = A;
            r.W = W;
            r.N = length(A);
            r.File = File;
            r.Path = Path;
            r.Labels = Labels;
            
        end
        function runREA(rea)

            rea.degree();
            rea.cluster();
            rea.weight();
            rea.eIndex();
            rea.laterality();
        end
        
        function [deg, mn, st] = degree(r)
    
            if isempty(r.Degree)
                A = r.A;
                r.Degree = sum(A~=0,2);
                r.MeanDegree = mean(r.Degree);
                r.StdDegree = std(r.Degree,1);
            end
            
            deg = r.Degree;
            mn = r.MeanDegree;
            st = r.StdDegree;
        end
        
        function [weiRea, weiIn, weiOut] = weight(r)
            
            if isempty(r.WeiDeg)
                A = r.A;
                W = r.W;

                r.WeiDeg = sum(A,2);
                r.WeiDegIn = sum(W,1)';
                r.WeiDegOut = sum(W,2);
            end
            
            weiRea = r.WeiDeg;
            weiIn = r.WeiDegIn;
            weiOut = r.WeiDegOut;
        end
        
        function cluster = cluster(r)
            
                A = r.A;
                numNodes = length(A);
                r.Cluster=zeros(numNodes,1);

                for node=1:numNodes
                    neigh=find(A(node,:));
                    KK=length(neigh);
                    if KK>=2
                        S=A(neigh,neigh);
                        r.Cluster(node)=sum(S(:)>0)/(KK*(KK-1));
                    end
                end
                cluster = r.Cluster;

        end
        
        function [Eindex, InsideCon, OutsideCon] = eIndex(r)
            
           
            if isempty(r.Eindex)
                 
                A = r.A;
                N = r.N;
                lab = r.Labels;
                
                
                H=zeros(N,2);
                
                if(length(lab)> 0 )
                
                    for i=1:N
                        SameSide=0;
                        OtherSide=0;
                        numSource=str2double(lab{i}(end));
                        if ~isnan(numSource)
                            nodeSourceSide=mod(numSource,2);
                            V=find(A(i,:));
                            KK=length(V);
                            
                            for j=1:KK
                                numTarget=str2double(lab{V(j)}(end));
                                if ~isnan(numTarget) %% if is a number
                                    nodeTargetSide=mod(numTarget,2);
                                
                                    if nodeTargetSide==nodeSourceSide
                                        SameSide=SameSide+A(i,V(j));
                                    else
                                        OtherSide=OtherSide+A(i,V(j));
                                    end
                                end
                            end
                        end
                        r.InsideCon(i)= SameSide ;
                        r.OutsideCon(i) = OtherSide; 

                    end
                    r.Eindex=(r.OutsideCon -r.InsideCon)./(r.InsideCon+r.OutsideCon);
                end
            end
            
            Eindex = r.Eindex;
            InsideCon =  r.InsideCon;
            OutsideCon = r.OutsideCon;
        end
        
        function [Lindex] = laterality(r)
            
           
            if isempty(r.Lindex)
                 
                W = r.W;
                N = r.N;
                lab = r.Labels;
                
                
                
                if(~isempty(lab))
                
                    for i=1:N
                        ll=0;
                        lr=0;
                        rl=0;
                        rr=0;
                        numSource=str2double(lab{i}(end));
                        if ~isnan(numSource)
                            nodeSourceSide=mod(numSource,2);
                            V=find(W(i,:));
                            KK=length(V);
                            
                            for j=1:KK
                                numTarget=str2double(lab{V(j)}(end));
                                if ~isnan(numTarget) %% if is a number
                                    nodeTargetSide=mod(numTarget,2);
                                    
                                    if nodeTargetSide ==1 && nodeSourceSide==1
                                        ll = ll + +W(i,V(j));
                                    elseif nodeTargetSide ==1 && nodeSourceSide==0
                                        lr = lr + +W(i,V(j));
                                    elseif nodeTargetSide ==0 && nodeSourceSide==1
                                        rl = rl + +W(i,V(j));
                                     elseif nodeTargetSide ==0 && nodeSourceSide==0
                                        rr = rr + +W(i,V(j));
                                    end
                                end
                            end
                        end
                        LL(i)=ll;
                        LR(i)=lr;
                        RL(i)=rl;
                        RR(i)=rr;
                      
                    end
                    r.Lindex= (((LL - RL)-(RR - LR))./(LL+LR+RL+RR));
                end
            end
            Lindex = r.Lindex;
           
        end
        
        function tabRea = reaCorrelationNoZeros(r)
            
            A = r.A;
            tabRea = table();

            reaIn = tril(A,-1);
            [target, source, weight] = find(reaIn);
            type = cell(length(target),1);
            type(:) = {'Undirected'};

            if ~isempty(source)
                tabRea = horzcat(num2cell(source),num2cell(target),num2cell(weight),type);
            end
        end
        
         function tabRea = reaCorrelation(r)
            
            A = r.A;
            tabRea = table();

            REATVG = A + 1;
            REA = triu(REATVG,1);
            [target, source, weight] = find(REA');
            weight = weight - 1;
            
            type = cell(length(target),1);
            type(:) = {'Undirected'};

            if ~isempty(source)      
                tabRea = horzcat(num2cell(source),num2cell(target),num2cell(weight),type);
            end
         end
        
         function tabRea = IRea(r)
             
             tabRea = horzcat(r.Labels,num2cell(r.Degree), num2cell(r.Cluster), num2cell(r.WeiDeg), num2cell(r.WeiDegIn), num2cell(r.WeiDegOut), num2cell(r.OutsideCon'), num2cell(r.InsideCon'), num2cell(r.Eindex'), num2cell(r.Lindex'));
         end
         
         function write (r)
            NameSync = 'MoS';
            REA_Tot = r.reaCorrelation(); %% with zero correlations
            REA_G = r.reaCorrelationNoZeros(); %% without zero correlations
            iREA = r.IRea();
            if ~isempty(r.Labels)
                idS=[REA_Tot{:,1}]';
                idT=[REA_Tot{:,2}]';
                REA_Tot(:,1)=r.Labels(idS);
                REA_Tot(:,2)=r.Labels(idT);
                if(~isempty(REA_G))
                    idS=[REA_G{:,1}]';
                    idT=[REA_G{:,2}]';
                    REA_G(:,1)=r.Labels(idS);
                    REA_G(:,2)=r.Labels(idT);
                end
            end
            
            
            writeFile( ...
                iREA, ...
                {'Label', 'DegreeREA', 'ClusREA', 'WeiDegREA',	'WeiDegIn',	'WeiDegOut',  'InsideCon', 'OutsideCon', 'Eindex', 'Lindex'}, ...
                [r.Path '\' r.File '_iREA_' NameSync], ...
                '.txt');
             writeFile( ...
                REA_Tot, ...
                {'Source', 'Target', 'Weight', 'Type'}, ...
                [r.Path '\' r.File '_REA_T_' NameSync], ...
                '.txt');
            if(~isempty(REA_G))
                writeFile( ...
                 REA_G, ...
                    {'Source', 'Target', 'Weight', 'Type'}, ...
                    [r.Path '\' r.File '_REA_G' NameSync], ...
                    '.txt');
            end

         end
        
    end
end