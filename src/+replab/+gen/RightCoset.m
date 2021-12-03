classdef RightCoset < replab.RightCoset & replab.gen.Coset

    methods

        function self = RightCoset(type, nice, niceIsomorphism, representative, subgroup, group)
            assert(isa(type, 'replab.gen.FiniteGroupType'));
            assert(isa(nice, 'replab.RightCoset'));
            assert(isa(niceIsomorphism, 'replab.gen.NiceIsomorphism'));
            self.type = type;
            self.nice = nice;
            self.niceIsomorphism = niceIsomorphism;
            self.representative_ = representative;
            self.subgroup = subgroup;
            self.group = group;
        end

    end

end
