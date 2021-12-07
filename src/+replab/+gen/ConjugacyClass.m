classdef ConjugacyClass < replab.ConjugacyClass & replab.gen.FiniteSet

    methods

        function self = ConjugacyClass(type, nice, niceIsomorphism, group)
            assert(isa(type, 'replab.gen.FiniteGroupType'));
            assert(isa(nice, 'replab.ConjugacyClass'));
            assert(isa(niceIsomorphism, 'replab.gen.NiceIsomorphism'));
            self.type = type;
            self.nice = nice;
            self.niceIsomorphism = niceIsomorphism;
            self.group = group;
        end

    end

    methods % Implementations

        % ConjugacyClass

        function o = elementOrder(self)
            o = self.nice.elementOrder;
        end

        function l = knownRepresentativeCentralizer(self)
            l = self.nice.knownRepresentativeCentralizer;
        end

        function G = representativeCentralizer(self)
            G = self.niceIsomorphism.preimageGroup(self.nice.representativeCentralizer);
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
