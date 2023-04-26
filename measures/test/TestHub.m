classdef TestHub < matlab.unittest.TestCase
    methods (Test)
        function testCalculateBinary(testCase)
            % Test the binary hub calculation
            adjMatrix = [0 1 1 0;
                         1 0 1 1;
                         1 1 0 1;
                         0 1 1 0];

            hub = Hub();
            result = hub.calculateBinary(adjMatrix);

            % The expected hub status vector for the given adjacency matrix
            expected = [0; 1; 1; 0];
            
            testCase.verifyEqual(result, expected);
        end

        function testCalculateWeighted(testCase)
            % Test the weighted hub calculation (it should return an empty array)
            adjMatrix = [0 2 3 0;
                         2 0 1 4;
                         3 1 0 2;
                         0 4 2 0];

            hub = Hub();
            result = hub.calculateWeighted(adjMatrix);

            % The expected result for weighted graphs is an empty array
            expected = [];
            
            testCase.verifyEqual(result, expected);
        end
    end
end