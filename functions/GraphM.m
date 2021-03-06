classdef GraphM < handle
   
    properties (Constant)
            
        DEGREE = 1;
        DEGREE_NAME = 'degree';
        DEGREE_NODAL = 'false';
        DEGREE_TXT = 'The degree of a node is the number of edges connceted to the node.';

        Name = { ...
           GraphM.DEGREE_NAME ...
        };

        NODAL = { ...
           GraphM.DEGREE_NODAL ...
        };
        TXT = { ...
            GraphM.DEGREE_TXT ...
        };
    end
    
    properties (GetAccess = public, SetAccess = protected)
        A % connection matrix
        P % coefficient p-values
        
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
        function g = GraphM(A)
            if nargin > 0
                g.A = double(A);
                g.N = length(g.A);
            end
            
        end
        function [deg, mn, st] = binDegree(g)
    
            if isempty(g.Degree)
                A = g.A~=0;
                g.Degree = sum(A,2);
                g.mn = mean(g.Degree);
                g.st = std(g.Degree,1);
            end
            
            deg = g.Degree;
            mn = g.mn;
            st = g.st;
        end
        
        function H = eIndex(g, lab)
            
            if isempty(g.H)
                A = g.A;
                N = g.N;
                
                g.H=zeros(N,2);
                
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
                                        SameSide=SameSide+1;
                                    else
                                        OtherSide=OtherSide+1;
                                    end
                                end
                            end
                        end
                        g.H(i,:)=[OtherSide,SameSide]; %% E-Index

                    end
                end
            end
            
            H = g.H;
        end
        
        function E = edges(g)
    
            if isempty(g.E)
                A = g.A;
                g.E = sum(sum(A))/2;
            end
            
            E = g.E;
        end
        
        function wei = weight(g)
    
            if isempty(g.Weight)
                A = g.A;
                A = A(:);
                A(A == 0) = [];
                g.Weight = A;
            end
            
            wei = g.Weight;
            
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
                g.dist = GraphM.removediagonal(g.dist);
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
        function r = diameter(g)
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
            
            r = max(g.eccentricity());
        end
        function [ecc,eccin,eccout] = eccentricity(g)
            % ECCENTRICITY eccentricity of nodes
            %
            % [ECC,ECCIN,ECCOUT] = ECCENTRICITY(G) calculates the eccentricity ECC,
            %   in-eccentricity ECCIN and out-eccentricity ECCOUT of all nodes in
            %   the graph G.
            %
            % The node eccentricity is the maximal shortest path length between a
            %   node and any other node. The eccentricity of a node is the maximum
            %   of the in and out-eccentricity of that node.
            %
            % Reference: "Combinatorics and Graph theory", J.M. Harris, J.L.Hirst
            %            and M.J. Mossinghoff
            %
            % See also Graph, distance, pl.
            
            if isempty(g.ecc) || isempty(g.eccin) || isempty(g.eccout)
                D = g.distance();
                D = GraphM.removediagonal(D,Inf);
                 
                g.eccin = max(D.*(D~=Inf),[],1);  % in-eccentricy = max of distance along column
                g.eccout = max(D.*(D~=Inf),[],2)';  % out-eccentricy = max of distance along column
                g.ecc = max(g.eccin,g.eccout);
            end
            
            ecc = g.ecc;
            eccin = g.eccin;
            eccout = g.eccout;
        end
        
        function [clo,cloin,cloout] = closeness(g)
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
        function [clNode] = cluster(g)
            
                A = g.A;
                numNodes = length(A);
                clNode=zeros(numNodes,1);

                for node=1:numNodes
                    neigh=find(A(node,:));
                    KK=length(neigh);
                    if KK>=2
                        S=A(neigh,neigh);
                        clNode(node)=sum(S(:))/(KK*(KK-1));
                    end
                end

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
        function [c,cin,cout] = pathLength(g)
            
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
            
            c = g.c;
            cin = g.cin;
            cout = g.cout;
        end
        
        function hubs = hub(g)
            
            if isempty(g.Degree)
                g.binDegree();
            end
            
            if isempty(g.hubs)
                A = g.A;
                
                g.hubs = g.Degree > (g.mn + 2* g.st);
            end
            
            hubs = g.hubs; 
            
        end
        
        function [GEdiff,Ediff] = diffusion_efficiency(g)
            % DIFFUSION_EFFICIENCY      Global mean and pair-wise diffusion efficiency
            %
            %   Withdrawn from Brain connectivity toolbox
            %
            %   [GEdiff,Ediff] = diffusion_efficiency(adj);
            %
            %   The diffusion efficiency between nodes i and j is the inverse of the
            %   mean first passage time from i to j, that is the expected number of
            %   steps it takes a random walker starting at node i to arrive for the
            %   first time at node j. Note that the mean first passage time is not a
            %   symmetric measure -- mfpt(i,j) may be different from mfpt(j,i) -- and
            %   the pair-wise diffusion efficiency matrix is hence also not symmetric.
            %
            %
            %   Input:
            %       adj,    Weighted/Unweighted, directed/undirected adjacency matrix
            %
            %
            %   Outputs:
            %       GEdiff, Mean Global diffusion efficiency (scalar)
            %       Ediff,  Pair-wise diffusion efficiency (matrix)
            %
            %
            %   References: Go??i J, et al (2013) PLoS ONE
            %
            %   Joaquin Go??i and Andrea Avena-Koenigsberger, IU Bloomington, 2012

            if isempty(g.GEdiff)
                A = g.A;
                display(A);
                n = size(A,1);
                mfpt = mean_first_passage_time(A);
                g.Ediff = 1./mfpt;
                g.Ediff(eye(n)>0) = 0;
                g.GEdiff = sum(g.Ediff(~eye(n)>0))/(n^2-n);
            end
            Ediff = g.Ediff;
            GEdiff = g.GEdiff;
        end
        
        function  [Erange,eta,Eshort,fs] = erange(g)
            %ERANGE     Shortcuts
            %
            %   Withdrawn from Brain connectivity toolbox
            %
            %   Shorcuts are central edges which significantly reduce the
            %   characteristic path length in the network.
            %
            %   Input:      CIJ,        binary directed connection matrix
            %
            %   Outputs:    Erange,     range for each edge, i.e. the length of the
            %                           shortest path from i to j for edge c(i,j) AFTER
            %                           the edge has been removed from the graph.
            %               eta         average range for entire graph.
            %               Eshort      entries are ones for shortcut edges.
            %               fs          fraction of shortcuts in the graph.
            %
            %   Follows the treatment of 'shortcuts' by Duncan Watts
            %
            %
            %   Olaf Sporns, Indiana University, 2002/2007/2008


            if isempty(g.Erange)
                A = g.A;
                N = size(A,1);
                K = length(nonzeros(A));
                g.Erange = zeros(N,N);
                [i,j] = find(A==1);

                for c=1:length(i)
                    CIJcut = A;
                    CIJcut(i(c),j(c)) = 0;
                    [~, D] = reachdist(CIJcut);
                    g.Erange(i(c),j(c)) = D(i(c),j(c));
                end;

                % average range (ignore Inf)
                g.eta = sum(g.Erange((g.Erange>0)&(g.Erange<Inf)))/length(g.Erange((g.Erange>0)&(g.Erange<Inf)));

                % Original entries of D are ones, thus entries of Erange
                % must be two or greater.
                % If Erange(i,j) > 2, then the edge is a shortcut.
                % 'fshort' is the fraction of shortcuts over the entire graph.

                g.Eshort = g.Erange>2;
                g.fs = length(nonzeros(g.Eshort))/K;
            end
            Erange = g.Erange;
            eta = g.eta;
            Eshort = g.Eshort;
            fs = g.fs;
            
        end
        
        function [Wq,twalk,wlq] = findwalks(g)
            %FINDWALKS      Network walks
            %
            %   Withdrawn from Brain connectivity toolbox
            %
            %   Walks are sequences of linked nodes, that may visit a single node more
            %   than once. This function finds the number of walks of a given length, 
            %   between any two nodes.
            %
            %   Input:      CIJ         binary (directed/undirected) connection matrix
            %
            %   Outputs:    Wq          3D matrix, Wq(i,j,q) is the number of walks
            %                           from 'i' to 'j' of length 'q'.
            %               twalk       total number of walks found
            %               wlq         walk length distribution as function of 'q'
            %
            %   Notes: Wq grows very quickly for larger N,K,q. Weights are discarded.
            %
            %   Algorithm: algebraic path count
            %
            %
            %   Olaf Sporns, Indiana University, 2002/2007/2008

            % ensure CIJ is binary...
            
            if isempty(g.Wq)
                A = g.A;
                A = double(A~=0);

                N = size(A,1);
                g.Wq = zeros(N,N,N);
                CIJpwr = A;
                g.Wq(:,:,1) = A;
                for q=2:N
                   CIJpwr = CIJpwr*A;
                   g.Wq(:,:,q) = CIJpwr;
                end;

                % total number of walks
                g.twalk = sum(sum(sum(g.Wq)));

                % walk length distribution
                g.wlq = reshape(sum(sum(g.Wq)),1,N);
            end
            
            Wq = g.Wq;
            twalk = g.twalk;
            wlq = g.wlq;
            
        end
        
        function [Ci,Q]=modularity_und(g,gamma)
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
        
        function [R,Nk,Ek] = rich_club_bd(g,varargin)
            %RICH_CLUB_BD        Rich club coefficients (binary directed graph)
            %
            %   R = rich_club_bd(CIJ)
            %   [R,Nk,Ek] = rich_club_bd(CIJ,klevel)
            %
            %   The rich club coefficient, R, at level k is the fraction of edges that
            %   connect nodes of degree k or higher out of the maximum number of edges
            %   that such nodes might share.
            %
            %   Input:      CIJ,        connection matrix, binary and directed
            %            klevel,        optional input argument. klevel sets the
            %                              maximum level at which the rich club
            %                              coefficient will be calculated. If klevel is
            %                              not included the the maximum level will be
            %                              set to the maximum degree of CIJ.
            %
            %   Output:       R,        vector of rich-club coefficients for levels
            %                              1 to klevel.
            %                Nk,        number of nodes with degree>k
            %                Ek,        number of edges remaining in subgraph with
            %                              degree>k
            %
            %   Reference: Colizza et al. (2006) Nat. Phys. 2:110.
            %
            %   Martijn van den Heuvel, University Medical Center Utrecht, 2011

            if isempty(g.R)
                A = g.A;
                N = size(A,1);                    %#ok<NASGU>

                % definition of "degree" as used for RC coefficients
                % degree is taken to be the sum of incoming and outgoing connectons
                if isempty(g.Degree)
                    [~,~,degree] = g.Degree();
                else
                    degree = g.Degree;
                end

                if nargin == 1
                    klevel = max(degree);
                elseif nargin == 2
                    klevel = varargin{1};
                elseif nargin > 2
                    error('number of inputs incorrect. Should be [CIJ], or [CIJ, klevel]')
                end

                g.R = zeros(1,klevel);
                g.Nk = zeros(1,klevel);
                g.Ek = zeros(1,klevel);
                for k = 1:klevel
                    SmallNodes=find(degree<=k);       %get 'small nodes' with degree <=k
                    subCIJ=A;                       %extract subnetwork of nodes >k by removing nodes <=k of CIJ
                    subCIJ(SmallNodes,:)=[];          %remove rows
                    subCIJ(:,SmallNodes)=[];          %remove columns
                    g.Nk(k)=size(subCIJ,2);             %number of nodes with degree >k
                    g.Ek(k)=sum(subCIJ(:));             %total number of connections in subgraph
                    g.R(k)=g.Ek(k)/(g.Nk(k)*(g.Nk(k)-1));     %unweighted rich-club coefficient
                end
            end
            
            R = g.R;
            Nk = g.Nk;
            Ek = g.Ek;
        end
    end
    
    
    methods (Static)
        function z = zeros(varargin)
            if(nargin == 0)
                z = GraphM;
            elseif any([varargin{:}] <= 0)
                z = GraphM.empty(varargin{:});
            else
                z = repmat(GraphM,varargin{:});
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