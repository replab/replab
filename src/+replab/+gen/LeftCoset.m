classdef LeftCoset < replab.LeftCoset & replab.gen.Coset

    methods

        function self = LeftCoset(type, nice, niceIsomorphism)
            self@replab.gen.FiniteSet(type, nice, niceIsomorphism);
        end

    end

end
