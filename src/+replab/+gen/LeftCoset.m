classdef LeftCoset < replab.LeftCoset & replab.gen.FiniteSet

    methods

        function self = LeftCoset(type, nice, niceIsomorphism)
            self@replab.gen.FiniteSet(type, nice, niceIsomorphism);
        end

    end

end
