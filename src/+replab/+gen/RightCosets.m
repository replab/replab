classdef RightCosets < replab.RightCosets

    properties (SetAccess = protected)
        nice % (`+replab.RightCosets`): Nice object where computations are done
        niceIsomorphism % (`+replab.gen.NiceIsomorphism`): Order-preserving isomorphism from a group containing ``self`` to a group containing `.nice`
    end

    methods

        function self = RightCosets(nice, niceIsomorphism, group, subgroup)
            assert(isa(nice, 'replab.RightCosets'));
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

        function T = transversal(self)
            T = cellfun(@(t) self.niceIsomorphism.preimageElement(t), self.nice.transversal, 'uniform', 0);
        end

    end

end
