classdef Collection < replab.Laws

    properties (SetAccess = protected)
        children % (cell(1,\*) of `+replab.Laws`): Laws instances in the collection
    end

    methods

        function self = Collection(children)
            self.children = children;
        end

        function n = nChildren(self)
            n = length(self.children);
        end

        function c = child(self, i)
            c = self.children{i};
        end

        function [testNames testFuns] = getTestCases(self)
            testNames = {};
            testFuns = {};
            for i = 1:self.nChildren
                c = self.child(i);
                [testNames1 testFuns1] = c.getTestCases;
                for j = 1:length(testNames1)
                    testNames1{j} = sprintf('%s(%d)', testNames1{j}, i);
                end
                testNames = horzcat(testNames, testNames1);
                testFuns = horzcat(testFuns, testFuns1);
            end
        end

    end

end
