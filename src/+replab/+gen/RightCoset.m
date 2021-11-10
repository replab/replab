classdef RightCoset < replab.RightCoset & replab.gen.FiniteSet

    methods

        function self = RightCoset(type, nice, niceIsomorphism)
            self@replab.gen.FiniteSet(type, nice, niceIsomorphism);
        end

    end

end
