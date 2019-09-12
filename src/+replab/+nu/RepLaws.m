classdef RepLaws < replab.Laws
    
    properties
        rep % non unitary representation being tested
        G % group of which rep is a representation
        L % Invertible matrices on R or C
        M % n x n matrices over R, or C
    end
    
    methods
        
        function self = RepLaws(rep)
            self.rep = rep;
            d = self.rep.dimension;
            self.G = rep.group;
            self.L = replab.domain.GeneralLinearMatrices(rep.field, d);
            switch rep.field
              case 'R'
                self.M = replab.domain.RealMatrices(d, d);
              case 'C'
                self.M = replab.domain.ComplexMatrices(d, d);
              otherwise
                error('Unknown field');
            end            
        end
        
        function morphismLaws = laws_asGroupHomomorphism(self)
            morphismLaws = replab.GroupMorphismLaws(@(g) self.rep.image(g), self.G, self.L);
        end
        
    end
end
