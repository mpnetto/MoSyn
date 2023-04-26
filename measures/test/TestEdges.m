classdef TestEdges < matlab.unittest.TestCase
    methods (Test)
        function testCalculateBinary(testCase)
            % Create an Edges object.
            edgesObj = Edges();
            
            % Define a binary adjacency matrix.
            adjMatrix = [0 1 1;
                         1 0 0;
                         1 0 0];
            
            % Expected value for the number of edges.
            expectedValue = 2;
            
            % Call the calculateBinary method.
            calculatedValue = edgesObj.calculateBinary(adjMatrix);
            
            % Verify that the expected and calculated values match.
            testCase.verifyEqual(calculatedValue, expectedValue);
        end
        
        function testCalculateWeighted(testCase)
            % Create an Edges object.
            edgesObj = Edges();
            
            % Define a weighted adjacency matrix.
            adjMatrix = [0 2 3;
                         2 0 0;
                         3 0 0];
            
            % Expected value for the number of edges.
            expectedValue = 2;
            
            % Call the calculateWeighted method.
            calculatedValue = edgesObj.calculateWeighted(adjMatrix);
            
            % Verify that the expected and calculated values match.
            testCase.verifyEqual(calculatedValue, expectedValue);
        end
    end
end