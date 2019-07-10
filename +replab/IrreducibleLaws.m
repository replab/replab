classdef IrreducibleLaws < replab.Laws
    
    properties
        irreducible; % Irreducible decomposition
        M;
    end
    
    methods
        
        function self = IrreducibleLaws(irreducible)
            self.irreducible = irreducible;
            d = self.irreducible.parent.dimension;
            switch irreducible.parent.field
              case 'R'
                self.M = replab.domain.RealMatrices(d, d);
              case 'C'
                self.M = replab.domain.ComplexMatrices(d, d);
              otherwise
                error('Unknown field');
            end
        end
        
        function law_decomposes_entire_space(self)
            A = self.irreducible.rep.A;
            self.M.assertEqv(A' * A, eye(self.irreducible.parent.dimension));
        end
        
        function isotypicLaws = laws_isotypic_components(self)
            children = cellfun(@(x) replab.IsotypicLaws(x), self.irreducible.components, 'uniform', 0);
            isotypicLaws = replab.LawsCollection(children);
        end
        
    end
end
