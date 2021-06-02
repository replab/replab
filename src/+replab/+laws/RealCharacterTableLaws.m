classdef RealCharacterTableLaws < replab.laws.CharacterTableLaws
% Law checks for real character tables

    methods

        function self = RealCharacterTableLaws(C)
            self@replab.laws.CharacterTableLaws(C);
        end

        function law_characters_are_real_(self)
            if self.C.group.order > 4096
                return
            end
            self.assert(isreal(self.C.values));
        end

        function law_irrep_images_are_real_(self)
            if self.C.group.order > 4096
                return
            end
            for i = 1:self.C.nIrreps
                if self.C.hasIrrep(i)
                    irrep = self.C.irrep(i);
                    for j = 1:self.C.group.nGenerators
                        img = irrep.image(self.C.group.generator(j), 'exact');
                        self.assert(isreal(img));
                    end
                end
            end
        end

    end

end
