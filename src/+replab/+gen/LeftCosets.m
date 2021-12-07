classdef LeftCosets < replab.LeftCosets

    properties (SetAccess = protected)
        nice % (`+replab.LeftCosets`): Nice object where computations are done
        niceIsomorphism % (`+replab.gen.NiceIsomorphism`): Order-preserving isomorphism from a group containing ``self`` to a group containing `.nice`
    end

    methods

        function self = LeftCosets(nice, niceIsomorphism, group, subgroup)
            assert(isa(nice, 'replab.LeftCosets'));
            self.nice = nice;
            self.niceIsomorphism = niceIsomorphism;
            self.group = group;
            self.subgroup = subgroup;
        end

    end

    methods % Implementations

        function t = cosetRepresentative(self, g)
            t = self.niceIsomorphism.preimageElement(self.nice.cosetRepresentative(self.niceIsomorphism.imageElement(g)));
        end

        function mu = leftAction(self)
            mu1 = self.nice.leftAction;
            mu = self.niceIsomorphism.restrictedSource(self.group).andThen(mu1);
        end

        function T = transversal(self)
            T = cellfun(@(t) self.niceIsomorphism.preimageElement(t), self.nice.transversal, 'uniform', 0);
        end

    end

end
