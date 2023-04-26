classdef TestClusterCoefficient < matlab.unittest.TestCase
    methods (Test)
        function testCalculateBinary(testCase)
            adjMatrix = [0 1 1;
                         1 0 1;
                         1 1 0];
            expected = [1; 1; 1];

            cc = ClusterCoefficient();
            result = cc.calculateBinary(adjMatrix);
            testCase.verifyEqual(result, expected);
        end
        
        function testCalculateWeighted(testCase)
            adjMatrix = [0 2 3;
                         2 0 2;
                         3 2 0];
            expected = [0.8; 1.5; 0.8];

            cc = ClusterCoefficient();
            result = cc.calculateWeighted(adjMatrix);
            testCase.verifyEqual(result, expected, 'AbsTol', 1e-6);
        end
    end
end