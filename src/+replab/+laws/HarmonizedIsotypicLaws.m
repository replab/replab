classdef HarmonizedIsotypicLaws < replab.laws.IsotypicLaws

    methods

        function self = HarmonizedIsotypicLaws(iso)
            self@replab.laws.IsotypicLaws(iso);
        end

        function law_all_irreps_same_basis_G(self, g)
            image1 = self.iso.irrep(1).image(g);
            for i = 2:self.iso.nIrreps
                imagei = self.iso.irrep(i).image(g);
                self.M1.assertEqv(image1, imagei);
            end
        end

    end

end
