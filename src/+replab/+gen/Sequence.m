classdef Sequence < replab.Sequence

    properties (SetAccess = protected)
        generic % (`+replab.Sequence`): Generic sequence
        genericIsomorphism % (`+replab.FiniteIsomorphism`): Isomorphism to the generic object
    end

    methods

        function self = Sequence(generic, genericIsomorphism)
            self@replab.Sequence(generic.nElements);
            self.generic = generic;
            self.genericIsomorphism = genericIsomorphism;

        end

    end

end
