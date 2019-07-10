classdef IsotypicLaws < replab.Laws
    
    properties
        iso; % Isotypic component
        C;   % commutant algebra
        G;   % group of which iso is a representation decomposition
        M1;   % n x n matrices over R or C
    end
    
    methods
        
        function self = IsotypicLaws(iso)
            assert(isa(iso, 'replab.Isotypic'));
            self.iso = iso;
            self.C = iso.rep.commutant;
            self.G = iso.parent.group;
            d1 = iso.copyDimension;
            switch iso.parent.field
              case 'R'
                self.M1 = replab.domain.RealMatrices(d1, d1);
              case 'C'
                self.M1 = replab.domain.ComplexMatrices(d1, d1);
              otherwise
                error('Unknown field');
            end
        end
        
        function irrepLaws = laws_subReps(self)
            children = cellfun(@(x) replab.IrrepLaws(x), self.iso.copies, 'uniform', 0);
            irrepLaws = replab.LawsCollection(children);
        end
        
        function law_all_irreps_equivalent_C(self, c)
            m = self.iso.multiplicity;
            d1 = self.iso.copyDimension;
            for i = 1:m
                for j = i:m
                    rows = (i-1)*d1+(1:d1);
                    cols = (j-1)*d1+(1:d1);
                    self.M1.assertNotEqv(c(rows, cols), zeros(d1, d1));
                end
            end
        end
        
        function law_all_irreps_same_basis_G(self, g)
            image1 = self.iso.copy(1).image(g);
            for i = 2:self.iso.nCopies
                imagei = self.iso.copy(i).image(g);
                self.M1.assertEqv(image1, imagei);
            end
        end
        
    end
end
