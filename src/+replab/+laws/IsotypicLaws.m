classdef IsotypicLaws < replab.Laws

    properties
        iso % replab.Isotypic: Isotypic component
        C % replab.Commutant: Commutant algebra
        G % replab.CompactGroup: group of which iso is a representation decomposition
        M1  % replab.domain.Matrices: n x n matrices over R or C
    end

    methods

        function self = IsotypicLaws(iso)
            assert(isa(iso, 'replab.Isotypic'));
            self.iso = iso;
            self.C = iso.commutant;
            self.G = iso.group;
            d1 = iso.irrepDimension;
            self.M1 = replab.domain.Matrices(iso.field, d1, d1);
        end

        function irrepLaws = laws_subReps(self)
            children = cellfun(@(x) x.laws, self.iso.irreps, 'uniform', 0);
            irrepLaws = replab.laws.Collection(children);
        end

        function law_all_irreps_equivalent_C(self, c)
            m = self.iso.multiplicity;
            d1 = self.iso.irrepDimension;
            for i = 1:m
                for j = i:m
                    rows = (i-1)*d1+(1:d1);
                    cols = (j-1)*d1+(1:d1);
                    self.M1.assertNotEqv(c(rows, cols), zeros(d1, d1));
                end
            end
        end

    end

end
