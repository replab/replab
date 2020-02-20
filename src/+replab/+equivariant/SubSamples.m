classdef SubSamples < replab.Str
% Describes a lazy linked list of samples derived from the samples of a parent equivariant space

    properties (SetAccess = protected)
        parent % (`+replab.+equivariant.Samples`): Parent representation samples
    end

    methods

        function self = SubSamples(E, parent)
            self@replab.equivariant.Samples(E);
            self.parent = parent;
        end

    end

    methods (Access = protected)

        function computeHead(self)
            [parentX parentErr] = self.parent.head;
            self.headSample_ = self.E.repR.F_internal * parentX * self.E.repC.H_internal;
            self.headError_ = NaN;
        end

        function computeTail(self)
            self.tail_ = replab.equivariant.SubSamples(self.E, self.parent.tail);
        end

    end

end
