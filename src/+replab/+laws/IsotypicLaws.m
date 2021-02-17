classdef IsotypicLaws < replab.laws.SubRepLaws

    methods

        function self = IsotypicLaws(iso)
            assert(isa(iso, 'replab.Isotypic'));
            self@replab.laws.SubRepLaws(iso);
        end

        function irrepLaws = laws_subReps(self)
            children = cellfun(@(x) x.laws, self.rep.irreps, 'uniform', 0);
            irrepLaws = replab.laws.Collection(children);
        end

        function law_all_irreps_same_basis_G(self, g)
            n = self.rep.nIrreps;
            for i = 1:n
                imagei = self.rep.irrep(i).image(g);
                for j = i+1:n
                    imagej = self.rep.irrep(j).image(g);
                    tol = self.rep.irrep(i).errorBound + self.rep.irrep(j).errorBound;
                    self.assertApproxEqual(imagei, imagej, tol);
                end
            end
        end

    end

end
