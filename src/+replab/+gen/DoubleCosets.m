classdef DoubleCosets < replab.DoubleCosets

    properties (SetAccess = protected)
        nice % (`+replab.DoubleCosets`): Nice object where computations are done
        niceIsomorphism % (`+replab.gen.NiceIsomorphism`): Order-preserving isomorphism from a group containing ``self`` to a group containing `.nice`
    end

    methods

        function self = DoubleCosets(nice, niceIsomorphism, group, leftSubgroup, rightSubgroup)
            assert(isa(nice, 'replab.DoubleCosets'));
            self.nice = nice;
            self.niceIsomorphism = niceIsomorphism;
            self.group = group;
            self.leftSubgroup = leftSubgroup;
            self.rightSubgroup = rightSubgroup;
        end

    end

    methods % Implementations

        function t = cosetRepresentative(self, g)
            t = self.niceIsomorphism.preimageElement(self.nice.cosetRepresentative(self.niceIsomorphism.imageElement(g)));
        end

        function T = transversal(self)
            T = cellfun(@(t) self.niceIsomorphism.preimageElement(t), self.nice.transversal, 'uniform', 0);
        end

    end

end
