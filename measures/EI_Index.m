
classdef EI_Index < Measure
    properties
        name = 'EI_Index';
        
        active = false;

        calculateType = 'Frame';

        supportedGraphs = ["GraphWD"];
    end
    
    methods
        function value = calculateBinary(~, adjMatrix)
            
       
        end
        function value = calculateWeighted(~, adjMatrix)           

            % internalConnections = 0;
            % externalConnections = 0;
            % 
            % N = length(adjMatrix);
            % 
            % for i = 1:N
            %     for j = 1:N
            %         if i ~= j && adjMatrix(i, j) > 0
            %             if ismember(i, group) && ismember(j, group)
            %                 internalConnections = internalConnections + adjMatrix(i, j);
            %             elseif ismember(i, group) || ismember(j, group)
            %                 externalConnections = externalConnections + adjMatrix(i, j);
            %             end
            %         end
            %     end
            % end
            % 
            % % Divide by 2 since we counted each internal edge twice
            % internalConnections = internalConnections / 2;
            % 
            % if internalConnections == 0
            %     ei_index = Inf;
            % else
            %     ei_index = externalConnections / internalConnections;
            % end
            value = [1,2,3];

        end
    end
end