classdef GroupSequence < replab.Sequence

    properties (SetAccess = protected)
        basis % (`+replab.+util.MixedRadix`): Basis to enumerate coset representatives
        lexChain % (`+replab.+bsgs.Chain`): Chain with a lexicographic base
    end

    methods

        function self = GroupSequence(nElements, lexChain)
            self@replab.Sequence(nElements);
            self.basis = replab.util.MixedRadix(lexChain.orbitSizes, 'oneBased', true, 'bigEndian', false);
            self.lexChain = lexChain;
        end

    end

    methods % Implementations

        % Sequence

        function el = at(self, ind)
            if ischar(ind)
                ind = vpi(ind);
            end
            el = self.lexChain.elementFromIndices(self.basis.ind2sub(ind));
        end

        function ind = find(self, el)
            ind = self.basis.sub2ind(self.lexChain.indicesFromElement(el));
        end

    end

end
