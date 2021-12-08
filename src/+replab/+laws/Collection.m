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

        function testCases = getTestCases(self, namePrefix, location)
            testCases = cell(1, 0);
            for i = 1:self.nChildren
                c = self.child(i);
                testCases = horzcat(testCases, c.getTestCases(sprintf('%s(%d)', namePrefix, i), location));
            end
        end

    end

end
