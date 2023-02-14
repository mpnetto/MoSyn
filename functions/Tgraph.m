classdef Tgraph < handle
   
    properties (Constant)
            
        DEGREE = 1;
        DEGREE_NAME = 'degree';
        DEGREE_NODAL = 'false';
        DEGREE_TXT = 'The degree of a node is the number of edges connceted to the node.';

        Name = { ...
           Tgraph.DEGREE_NAME ...
        };

        NODAL = { ...
           Tgraph.DEGREE_NODAL ...
        };
        TXT = { ...
            Tgraph.DEGREE_TXT ...
        };
    end
    
    properties (GetAccess = public, SetAccess = protected)
        A % connection matrix
        P % coefficient p-values

        time
        nodesLocation
        
        %N = nodeNumber;
        N
        
        % wei = weight(g)
        Weight
        
        %D = distance(g);
        D
        
        % [degree, mean] = binDegree(A)
        Degree
        mn
        st
        
        % [ecc,eccin,eccout] = eccentricity(g)
        ecc
        eccin
        eccout
        
        % [clo,cloin,cloout] = closeness(g)
        clo
        cloin
        cloout
        
        % b = betweenness(g,normalized)
        b
        
        % [cl, clNode] = cluster(g)
        clNode
        
        % dist = distance(A)
        dist
        
        % [c,cin,cout] = pathLength(g)
        c
        cin
        cout
        
        % hubs = hub(g)
        hubs
        
        % E = edges(g)
        E
        
        % H = eIndex(g, lab)
        H
        
        % [GEdiff,Ediff] = diffusion_efficiency(g)
        GEdiff
        Ediff
        
        % [Erange,eta,Eshort,fs] = erange(g)
        Erange
        eta
        Eshort
        fs
        
        % [Wq,twalk,wlq] = findwalks(CIJ)
        Wq
        twalk
        wlq
        
        % [Ci,Q]=modularity_und(A,gamma)
        Ci
        Q
        
        % [R,Nk,Ek] = rich_club_bd(CIJ,varargin)
        R
        Nk
        Ek
        
    end
   
    methods 
        
        function g = Tgraph(A, nodeLocation, time)
            if nargin > 0
                g.A = double(A);
                g.N = length(g.A);
            end

            if nargin > 1
                g.nodesLocation = nodeLocation;
                g.time = time;
            end
            
        end
        
        function d = degree(g)
    
            if isempty(g.Degree)
                A = g.A~=0;
                g.Degree = sum(A,2);
                g.mn = mean(g.Degree);
                g.st = std(g.Degree,1);
            end
            
            d = g.Degree;
            mn = g.mn;
            st = g.st;
        end
        
        function p = pathLength(g)
            
            if isempty(g.c) 
                A = g.A;
    
                if isempty(g.dist)
                    g.distance();
                end
                dist = g.dist;
                numNodes = length(A);

                g.cin = zeros(1,numNodes);
                for node = 1:1:numNodes
                    distNode = dist(:,node);
                    g.cin(node) = sum(distNode(distNode~=Inf))/length(nonzeros(distNode~=Inf));
                end

                g.cout = zeros(1,numNodes);
                for node = 1:1:numNodes
                    distNode = dist(node,:);
                    g.cout(node) = sum(distNode(distNode~=Inf))/length(nonzeros(distNode~=Inf));
                end

                g.c = mean([g.cin; g.cout],1);
            end
            
            p = g.c;
            cin = g.cin;
            cout = g.cout;
        end
        
        function dist = distance(g)
            
            if isempty(g.dist)
                A = g.A;
                
                len = 1;  
                g.dist = g.A;

                Lpath = g.A;
                Idx = true;
                while any(Idx(:))
                    len = len+1;
                    Lpath = Lpath*g.A;
                    Idx = (Lpath~=0)&(g.dist==0);
                    g.dist(Idx) = len;
                end

                g.dist(~g.dist) = inf;
                g.dist = Tgraph.removediagonal(g.dist);
            end
            dist = g.dist;
        end
        
        function r = radius(g)
            % RADIUS radius of a graph
            %
            % R = RADIUS(G) calculates the radius R of the graph G.
            %
            % Radius is calculated as the minimum eccentricity.
            %
            % Reference: "Combinatorics and Graph theory", J.M. Harris, J.L.Hirst
            %            and M.J. Mossinghoff
            %
            % See also Graph, diameter, eccentricity.
            
            r = min(g.eccentricity());
        end
        
        function d = diameter(g)
            % DIAMETER diameter of a graph
            %
            % R = DIAMETER(G) calculates the diameter R of the graph G.
            %
            % Diameter is calculated as the maximum eccentricity.
            %
            % Reference: "Combinatorics and Graph theory", J.M. Harris, J.L.Hirst
            %            and M.J. Mossinghoff
            %
            % See also Graph, radius, eccentricity.
            
            d = max(g.eccentricity());
        end
    
        function clu = cluster(g)
            
            A = g.A;
            numNodes = length(A);
            clNode=zeros(numNodes,1);

            for node=1:numNodes
                neigh=find(A(node,:));
                KK=length(neigh);
                if KK>=2
                    S=A(neigh,neigh);
                    clu(node)=sum(S(:))/(KK*(KK-1));
                end
            end

        end
        
        function Q = modularity(g,gamma)
            %MODULARITY_UND     Optimal community structure and modularity
            %
            %   Withdrawn from Brain connectivity toolbox
            %
            %   [Ci Q] = modularity_und(W,gamma);
            %
            %   The optimal community structure is a subdivision of the network into
            %   nonoverlapping groups of nodes in a way that maximizes the number of
            %   within-group edges, and minimizes the number of between-group edges.
            %   The modularity is a statistic that quantifies the degree to which the
            %   network may be subdivided into such clearly delineated groups.
            %
            %   Inputs:
            %       W,
            %           undirected weighted/binary connection matrix
            %       gamma,
            %           resolution parameter (optional)
            %               gamma>1,        detects smaller modules
            %               0<=gamma<1,     detects larger modules
            %               gamma=1,        classic modularity (default)
            %
            %   Outputs:    
            %       Ci,     optimal community structure
            %       Q,      maximized modularity
            %
            %   Note:
            %       This algorithm is essentially deterministic. The only potential
            %       source of stochasticity occurs at the iterative finetuning step, in
            %       the presence of non-unique optimal swaps. However, the present
            %       implementation always makes the first available optimal swap and
            %       is therefore deterministic.
            %
            %   References: 
            %       Newman (2006) -- Phys Rev E 74:036104, PNAS 23:8577-8582.
            %       Reichardt and Bornholdt (2006) Phys Rev E 74:016110.
            %
            %   2008-2016
            %   Mika Rubinov, UNSW
            %   Jonathan Power, WUSTL
            %   Dani Bassett, UCSB
            %   Xindi Wang, Beijing Normal University
            %   Roan LaPlante, Martinos Center for Biomedical Imaging

            %   Modification History:
            %   Jul 2008: Original (Mika Rubinov)
            %   Oct 2008: Positive eigenvalues made insufficient for division (Jonathan Power)
            %   Dec 2008: Fine-tuning made consistent with Newman's description (Jonathan Power)
            %   Dec 2008: Fine-tuning vectorized (Mika Rubinov)
            %   Sep 2010: Node identities permuted (Dani Bassett)
            %   Dec 2013: Gamma resolution parameter included (Mika Rubinov)
            %   Dec 2013: Detection of maximum real part of eigenvalues enforced (Mika Rubinov)
            %               Thanks to Mason Porter and Jack Setford, University of Oxford
            %   Dec 2015: Single moves during fine-tuning enforced (Xindi Wang)
            %   Jan 2017: Removed node permutation and updated documentation (Roan LaPlante)

            if isempty(g.Ci)
                if ~exist('gamma','var')
                    gamma = 1;
                end
                
                A = g.A;
                N=length(A);                            %number of vertices
                % n_perm = randperm(N);                   %DB: randomly permute order of nodes
                % A = A(n_perm,n_perm);                   %DB: use permuted matrix for subsequent analysis
                K=sum(A);                               %degree
                m=sum(K);                               %number of edges (each undirected edge is counted twice)
                B=A-gamma*(K.'*K)/m;                    %modularity matrix
                g.Ci=ones(N,1);                           %community indices
                cn=1;                                   %number of communities
                U=[1 0];                                %array of unexamined communites

                ind=1:N;
                Bg=B;
                Ng=N;

                while U(1)                              %examine community U(1)
                    Bg(isnan(Bg))=0;
                    Bg(isinf(Bg))=0;
                    [V,D]=eig(Bg);
                    [~,i1]=max(real(diag(D)));          %maximal positive (real part of) eigenvalue of Bg
                    v1=V(:,i1);                         %corresponding eigenvector

                    S=ones(Ng,1);
                    S(v1<0)=-1;
                    q=S.'*Bg*S;                         %contribution to modularity

                    if q>1e-10                       	%contribution positive: U(1) is divisible
                        qmax=q;                         %maximal contribution to modularity
                        Bg(logical(eye(Ng)))=0;      	%Bg is modified, to enable fine-tuning
                        indg=ones(Ng,1);                %array of unmoved indices
                        Sit=S;
                        while any(indg)                 %iterative fine-tuning
                            Qit=qmax-4*Sit.*(Bg*Sit); 	%this line is equivalent to:
                            [qmax,imax]=max(Qit.*indg); %for i=1:Ng
                            Sit(imax)=-Sit(imax);       %	Sit(i)=-Sit(i);
                            indg(imax)=nan;             %	Qit(i)=Sit.'*Bg*Sit;
                            if qmax>q                   %	Sit(i)=-Sit(i);
                                q=qmax;                 %end
                                S=Sit;
                            end
                        end

                        if abs(sum(S))==Ng              %unsuccessful splitting of U(1)
                            U(1)=[];
                        else
                            cn=cn+1;
                            g.Ci(ind(S==1))=U(1);         %split old U(1) into new U(1) and into cn
                            g.Ci(ind(S==-1))=cn;
                            U=[cn U];                   %#ok<AGROW>
                        end
                    else                                %contribution nonpositive: U(1) is indivisible
                        U(1)=[];
                    end

                    ind=find(g.Ci==U(1));                 %indices of unexamined community U(1)
                    bg=B(ind,ind);
                    Bg=bg-diag(sum(bg));                %modularity matrix for U(1)
                    Ng=length(ind);                     %number of vertices in U(1)
                end

                s=g.Ci(:,ones(1,N));                      %compute modularity
                g.Q=~(s-s.').*B/m;
                g.Q=sum(g.Q(:));
            end
            Ci=g.Ci;
            Q=g.Q;
            
        end
       
        function clo = closeness(g)
            % CLOSENESS closeness centrality of nodes
            %
            % [CLO,CLOIN,CLOOUT] = CLOSENESS(G) calculates the closeness centrality CLO,
            %   in-closeness centrality CLOIN and out-closeness centrality CLOOUT of all
            %   nodes in the graph G.
            %
            % The closeness centrality is the inverse of the average shortest path lengths
            %   of one node to all other nodes.
            %
            % See also Graph, pl.
            
            if isempty(g.clo)
                [c,cin,cout] = g.pathLength();
                g.clo = c.^-1;
                g.cloin = cin.^-1;
                g.cloout = cout.^-1;
            end

            clo = g.clo;
            cloin = g.cloin;
            cloout = g.cloout;
        end

        function b = betweenness(g,normalized)
            % BETWEENNESS betweenness centrality of a node
            %
            % B = BETWEENNESS(G) calculates the betweenness centrality B of all
            %   nodes of the graph G.
            %
            % B = BETWEENNESS(G,NORMALIZED) calculates the betweenness centrality
            %   B and normalizes it to the range [0,1] if NORMALIZED is equal to true (default).
            %   If NORMALIZED is equal to false, the betweenness is not normalized.
            %
            % Betweenness centrality of a node is defined as the fraction of
            %   all shortest paths in a graph G that contain that node. Nodes
            %   that have high value of betweenness centrality participate in
            %   a large number of shortest paths.
            %
            % Betweenness centrality B may be normalised to the range [0,1] as
            %   B/[(N-1)(N-2)], where N is the number of nodes in the graph.
            %
            %   Reference: "Betweenness centrality: Algorithms and Lower Bounds", S.Kintali
            %              (generalization to directed and disconnected graphs)
            %              "Complex network measures of brain connectivity: Uses and
            %              interpretations", M.Rubinov and O.Sporns
            %
            % See also GraphBD.
            
            if isempty(g.b)
                A = g.A;
                N = g.N;  % number of nodes
                I = eye(N)~=0;  % logical identity matrix
                d = 1;  % path length
                NPd = A;  % number of paths of length |d|
                NSPd = NPd;  % number of shortest paths of length |d|
                NSP = NSPd; NSP(I) = 1;  % number of shortest paths of any length
                L = NSPd; L(I) = 1;  % length of shortest paths
                
                % calculate NSP and L
                while find(NSPd,1);
                    d = d+1;
                    NPd = NPd*A;
                    NSPd = NPd.*(L==0);
                    NSP = NSP+NSPd;
                    L = L+d.*(NSPd~=0);
                end
                L(~L) = inf; L(I) = 0;  % L for disconnected vertices is inf
                NSP(~NSP) = 1;  % NSP for disconnected vertices is 1
                
                At = A.';
                DP = zeros(N);  % vertex on vertex dependency
                diam = d-1;  % graph diameter
                % calculate DP
                for d = diam:-1:2
                    DPd1 = (((L==d).*(1+DP)./NSP)*At).*((L==(d-1)).*NSP);
                    DP = DP + DPd1;  % DPd1: dependencies on vertices |d-1| from source
                end
                
                g.b = sum(DP,1);  % compute betweenness
                
            end
            b = g.b;
            
            if exist('normalized') && normalized
                N = g.nodenumber();
                b = b/((N-1)*(N-2));
            end
        end
        
 
    end
   
    methods (Static)
        function z = zeros(varargin)
            if(nargin == 0)
                z = Tgraph;
            elseif any([varargin{:}] <= 0)
                z = Tgraph.empty(varargin{:});
            else
                z = repmat(Tgraph,varargin{:});
            end
        end
        function B = removediagonal(A, value)
            % REMOVEDIAGONAL replaces matrix diagonal with given value
            %
            % B = REMOVEDIAGONAL(A,VALUE) sets the values in the diagonal of any
            %   matrix A to VALUE and returns the resulting matrix B.
            %   If VALUE is not inputted, the default is 0.
            %
            % See also Graph.
            
            if nargin < 2
                value = 0;
            end
            
            B = A;
            B(1:1+length(A):end) = value; 
        end
    end

end