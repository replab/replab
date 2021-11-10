classdef FiniteSet < replab.FiniteSet
% Describes a subset of a finite group using an isomorphism into a nice group where computations are easier
%
% We require that ``self.niceIsomorphism.source`` contains ``self``, and that ``self.niceIsomorphism`` preserves the total order
% of ``self.type`` and ``self.nice.type``.

    properties (SetAccess = protected)
        nice % (`+replab.FiniteSet`): Nice object where computations are done
        niceIsomorphism % (`+replab.gen.NiceIsomorphism`): Order-preserving isomorphism from a group containing ``self`` to a group containing `.nice`
    end

    methods

        function self = FiniteSet(type, nice, niceIsomorphism)
            assert(isa(type, 'replab.gen.FiniteGroupType'));
            assert(isa(nice, 'replab.FiniteSet'));
            assert(isa(niceIsomorphism, 'replab.gen.NiceIsomorphism'));
            self.type = type;
            self.nice = nice;
            self.niceIsomorphism = niceIsomorphism;
        end

    end

    methods % Implementations

        function b = contains(self, el)
            b = self.nice.contains(self.niceIsomorphism.imageElement(el));
        end

        function E = elements(self)
            E = cellfun(@(e) self.niceIsomorphism.preimageElement(e), self.nice.elements, 'uniform', 0);
        end

        function E = elementsSequence(self)
            E = self.nice.elementsSequence.imap(self.niceIsomorphism.inverse);
        end

        function s = nElements(self)
            s = self.nice.nElements;
        end

        function r = representative(self)
            self.representative_ = niceIsomorphism.preimageElement(nice.representative);
        end

        function s = setProduct(self)
            s = self.nice.setProduct.imap(self.niceIsomorphism.inverse);
        end

    end

end
