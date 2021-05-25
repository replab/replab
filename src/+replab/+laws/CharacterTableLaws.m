classdef CharacterTableLaws < replab.Laws
% Law checks for character tables

    properties (SetAccess = protected)
        C % (`+replab.CharacterTable`): Character table
    end

    methods

        function self = CharacterTableLaws(C)
            self.C = C;
        end

        function laws = laws_classes(self)
            laws = replab.laws.Collection({self.C.classes.laws});
        end

        function law_irrep_traces_match_character_values_(self)
            if self.C.group.order > 4096
                return
            end
            for i = 1:self.C.nIrreps
                irrep = self.C.irreps{i};
                if ~isempty(irrep)
                    for j = 1:self.C.nClasses
                        [c1, err1] = doubleApproximation(self.C.values(i, j));
                        c2 = trace(irrep.image(self.C.classes.classes{j}.representative));
                        self.assertApproxEqual(c1, c2, err1 + irrep.errorBound);
                    end
                end
            end
        end

        function law_orthogonality_rows_(self)
            n = self.C.nIrreps;
            res = replab.cyclotomic.zeros(n, n);
            for i1 = 1:n
                chi1 = self.C.character(i1);
                for i2 = 1:n
                    chi2 = self.C.character(i2);
                    res(i1, i2) = dot(chi1, chi2);
                end
            end
            self.assert(all(all(res == replab.cyclotomic.eye(n))));
        end

    end

end
