classdef TestPathLength < matlab.unittest.TestCase
    methods (Test)
        
        % Test case for binary graph.
        function testBinaryGraph(testCase)
            % Define the adjacency matrix for a binary graph.
            adjMatrix = [0 1 0 1;
                         1 0 1 0;
                         0 1 0 1;
                         1 0 1 0];
            
            % Create an instance of the PathLength class and calculate the path lengths.
            pl = PathLength();
            values = pl.calculateBinary(adjMatrix);
            
            % Define the expected path length values.
            expectedValues = [1.6667; 1.6667; 1.6667; 1.6667];
            
            % Verify that the calculated path lengths match the expected values.
            testCase.verifyEqual(values, expectedValues, 'RelTol', 1e-4);
        end

        % Test case for weighted graph.
        function testWeightedGraph(testCase)
            % Define the adjacency matrix for a weighted graph.
            adjMatrix = [0 2 0 3;
                         2 0 1 0;
                         0 1 0 4;
                         3 0 4 0];
            
            % Create an instance of the PathLength class and calculate the path lengths.
            pl = PathLength();
            values = pl.calculateWeighted(adjMatrix);
            
            % Define the expected path length values.
            expectedValues = [1.3333; 1.3333; 1.6667; 1.6667];
            
            % Verify that the calculated path lengths match the expected values.
            testCase.verifyEqual(values, expectedValues, 'RelTol', 1e-4);
        end
    end
end