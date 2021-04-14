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

        function law_orthogonality_columns_(self)
            if self.C.overC
                n = self.C.nClasses;
                for i1 = 1:n
                    for i2 = 1:n
                        r = dot(self.C.values(:, i1), self.C.values(:, i2));
                        if i1 ~= i2
                            self.assert(r == 0);
                        else
                            self.assert(r == replab.cyclotomic(self.C.classes.classes{i1}.representativeCentralizer.order));
                        end
                    end
                end
            end
        end

        function law_group_order_(self)
            if self.C.overC
                col = self.C.values(:, self.C.identityConjugacyClassIndex);
                order1 = sum(col.*col);
                order2 = replab.cyclotomic(self.C.group.order);
                self.assert(order1 == order2);
            end
        end

        function law_commutator_subgroup_(self)
        % The commutator subgroup of the group is the intersection of the kernels of the linear characters
            if self.C.overC
                chars = arrayfun(@(i) self.C.character(i), self.C.linearCharacterIndices, 'uniform', 0);
                kernels = cellfun(@(c) c.kernel, chars, 'uniform', 0);
                K = kernels{1};
                for i = 2:length(kernels)
                    K = K.intersection(kernels{i});
                end
                self.assert(K == self.C.group.derivedSubgroup);
            end
        end

    end

end
