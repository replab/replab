classdef LeftCoset < replab.LeftCoset & replab.gen.FiniteSet

    methods

        function self = LeftCoset(type, nice, niceIsomorphism, representative, subgroup, group)
            assert(isa(type, 'replab.gen.FiniteGroupType'));
            assert(isa(nice, 'replab.LeftCoset'));
            assert(isa(niceIsomorphism, 'replab.gen.NiceIsomorphism'));
            self.type = type;
            self.nice = nice;
            self.niceIsomorphism = niceIsomorphism;
            self.representative_ = representative;
            self.subgroup = subgroup;
            self.group = group;
        end

    end

    methods % Implementations

        % Str

        function names = hiddenFields(self)
            names = hiddenFields@replab.gen.FiniteSet(self);
        end

        % Obj

        function l = laws(self)
            l = laws@replab.LeftCoset(self);
        end

        % Domain

        function s = sample(self)
            s = sample@replab.LeftCoset(self);
        end

        % FiniteSet

        function C = imap(self, f)
            C = imap@replab.LeftCoset(self, f);
        end

        % gen.FiniteSet

        function l = compatibleWithNiceIsomorphism(self, iso)
            l = self.group.compatibleWithNiceIsomorphism(iso) && self.subgroup.compatibleWithNiceIsomorphism(iso);
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
