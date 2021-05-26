classdef ComplexCharacterTableLaws < replab.laws.CharacterTableLaws
% Law checks for complex character tables

    methods

        function self = ComplexCharacterTableLaws(C)
            self@replab.laws.CharacterTableLaws(C);
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
