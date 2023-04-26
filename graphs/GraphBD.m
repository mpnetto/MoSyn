classdef  GraphBD < Graph
   
    properties 
    end
    
    methods
        function obj = GraphBD(adjMatrix)
            if nargin > 0
                obj.AdjacencyMatrix = adjMatrix;
            end
        end
       
    end

    methods (Static)
        function z = zeros(varargin)
            if(nargin == 0)
                z = Subject;
            elseif any([varargin{:}] <= 0)
                z = GraphBD.empty(varargin{:});
            else
                z = repmat(GraphBD,varargin{:});
            end
        end
    end

end