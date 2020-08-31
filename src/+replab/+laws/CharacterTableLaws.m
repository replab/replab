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

        function law_irrep_(self)
            for i = 1:self.C.nIrreps
                irrep = self.C.irreps{i};
                if ~isempty(irrep)
                    charI = replab.Character.fromApproximateRep(irrep);
                    self.assert(charI == self.C.character(i));
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
            n = self.C.nClasses;
            for i1 = 1:n
                for i2 = 1:n
                    r = dot(self.C.characters(:, i1), self.C.characters(:, i2));
                    if i1 ~= i2
                        self.assert(r == 0);
                    else
                        self.assert(r == replab.cyclotomic.fromVPIs(self.C.classes.classes{i1}.representativeCentralizer.order));
                    end
                end
            end
        end

        function law_group_order_(self)
            col = self.C.characters(:, self.C.identityConjugacyClassIndex);
            order1 = sum(col.*col);
            order2 = replab.cyclotomic.fromVPIs(self.C.group.order);
            self.assert(order1 == order2);
        end

        function law_commutator_subgroup_(self)
        % The commutator subgroup of the group is the intersection of the kernels of the linear characters
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