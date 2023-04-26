classdef TestWeight < matlab.unittest.TestCase
    methods (Test)
        function testCalculateWeighted(testCase)
            % Create an Edges object.
            weightObj = Weight();
            
            % Define a binary adjacency matrix.
            adjMatrix = [0 2 3;
                         2 0 1;
                         3 1 0];
            
            % Calculate the sum of weights for each node using the Weight class.
            result = weightObj.calculateWeighted(adjMatrix);
            
            % Define the expected result.
            expected = [5; 3; 4];
            
            % Verify if the calculated result matches the expected result.
            testCase.verifyEqual(result, expected);
        end
        
    end
end