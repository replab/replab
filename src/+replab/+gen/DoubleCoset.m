classdef DoubleCoset < replab.DoubleCoset & replab.gen.FiniteSet

    methods

        function self = DoubleCoset(type, nice, niceIsomorphism)
            self@replab.gen.FiniteSet(type, nice, niceIsomorphism);
        end

    end

    methods % Implementations

        function l = compatibleWithNiceIsomorphism(self, iso)
            l = self.group.compatibleWithNiceIsomorphism(iso) && ...
                self.leftSubgroup.compatibleWithNiceIsomorphism(iso) && ...
                self.rightSubgroup.compatibleWithNiceIsomorphism(iso);
        end

    end

end
