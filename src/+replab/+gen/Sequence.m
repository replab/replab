classdef Sequence < replab.Sequence

    properties (SetAccess = protected)
        nice % (`+replab.Sequence`): Generic sequence
        niceIsomorphism % (`+replab.gen.NiceIsomorphism`): Isomorphism to the generic object
    end

    methods

        function self = Sequence(nice, niceIsomorphism)
            self@replab.Sequence(nice.nElements);
            self.nice = nice;
            self.niceIsomorphism = niceIsomorphism;
        end

        function obj = at(self, ind)
            obj = self.niceIsomorphism.preimageElement(self.nice.at(ind))
        end

        function ind = find(self, obj)
            if ~self.niceIsomorphism.sourceContains(obj)
                ind = vpi(0);
                return
            end
            ind = self.nice.find(self.niceIsomorphism.imageElement(obj));
        end

    end

end
