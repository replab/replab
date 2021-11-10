classdef GroupSequence < replab.Sequence

    properties (SetAccess = protected)
        basis % (`+replab.+util.MixedRadix`): Basis to enumerate coset representatives
        lexChain % (`+replab.+bsgs.Chain`): Chain with a lexicographic base
    end

    methods

        function self = GroupSequence(nElements, basis, lexChain)
            self@replab.Sequence(nElements);
            self.basis = basis;
            self.lexChain = lexChain;
        end

    end

    methods % Implementations

        % Sequence

        function el = at(self, ind)
            el = self.lexChain.elementFromIndices(self.basis.ind2sub(ind));
        end

        function ind = find(self, el)
            ind = basis.sub2ind(self.lexChain.indicesFromElement(el));
        end

    end

end
