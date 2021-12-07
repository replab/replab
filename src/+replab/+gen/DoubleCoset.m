classdef DoubleCoset < replab.DoubleCoset & replab.gen.FiniteSet

    methods

        function self = DoubleCoset(type, nice, niceIsomorphism, representative, leftSubgroup, rightSubgroup, group)
            assert(isa(type, 'replab.gen.FiniteGroupType'));
            assert(isa(nice, 'replab.DoubleCoset'));
            assert(isa(niceIsomorphism, 'replab.gen.NiceIsomorphism'));
            self.type = type;
            self.nice = nice;
            self.niceIsomorphism = niceIsomorphism;
            self.representative_ = representative;
            self.leftSubgroup = leftSubgroup;
            self.rightSubgroup = rightSubgroup;
            self.group = group;
        end

    end

    methods % Implementations

        function l = compatibleWithNiceIsomorphism(self, iso)
            l = self.group.compatibleWithNiceIsomorphism(iso) && ...
                self.leftSubgroup.compatibleWithNiceIsomorphism(iso) && ...
                self.rightSubgroup.compatibleWithNiceIsomorphism(iso);
        end

    end

    methods

        % Bugfix for Octave method selection

        function b = contains(self, el)
            b = contains@replab.gen.FiniteSet(self, el);
        end

        function E = elements(self)
            E = elements@replab.gen.FiniteSet(self);
        end

        function E = elementsSequence(self)
            E = elementsSequence@replab.gen.FiniteSet(self);
        end

        function s = nElements(self)
            s = nElements@replab.gen.FiniteSet(self);
        end

        function r = representative(self)
            r = representative@replab.gen.FiniteSet(self);
        end

        function S = setProduct(self)
            S = setProduct@replab.gen.FiniteSet(self);
        end

    end

end
