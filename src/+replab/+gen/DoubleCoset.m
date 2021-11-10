classdef DoubleCoset < replab.DoubleCoset & replab.gen.FiniteSet

    methods

        function self = DoubleCoset(type, nice, niceIsomorphism)
            self@replab.gen.FiniteSet(type, nice, niceIsomorphism);
        end

    end

end
