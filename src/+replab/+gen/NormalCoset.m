classdef NormalCoset < replab.NormalCoset & replab.gen.FiniteSet

    methods

        function self = NormalCoset(type, nice, niceIsomorphism)
            self@replab.gen.FiniteSet(type, nice, niceIsomorphism);
        end

    end

end
