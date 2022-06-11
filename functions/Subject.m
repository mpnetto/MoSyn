classdef Subject < handle & matlab.mixin.Copyable
   
    properties (Constant)
        
        DEGREE = 1;
        DEGREE_NAME = 'degree';
        DEGREE_NODAL = 'false';
        DEGREE_TXT = 'The degree of a node is the number of edges connceted to the node.';
       
        Name = { ...
           Subject.DEGREE_NAME ...
        };

        NODAL = { ...
           Subject.DEGREE_NODAL ...
        };
        TXT = { ...
            Subject.DEGREE_TXT ...
        };
    end
    
    properties (GetAccess = public, SetAccess = protected)
        A % connection matrix
        P % coefficient p-values
    end
    
    properties (Access = protected)
        
        %N = nodeNumber;
        N
        
        %D = distance(g);
        D
        
        % [deg, indeg, outdeg] = degree(g)
        deg
        indeg
        outdeg
    end
    
    
    methods (Access = protected)
        function g = GraphM(A)
            g.A = A;
            g.N = length(g.A);
        end
        function [degree, mean, std] = binDegree(A)

            degree = sum(A,1);
            mean = mean(degree);
            std = std(degree);
            
        end
        
    end
   
  
        
    
end