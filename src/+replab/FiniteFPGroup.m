classdef FiniteFPGroup < replab.FPGroup & replab.NiceFiniteGroup
% Describes a finite group with a presentation

    methods

        function setPermutationImages(self, permutations)
        % Sets the permutation realization of this group
        %
        % Enables to skip the Todd-Coxeter procedure when operating over the group.
            self.niceGroup_ = replab.PermutationGroup.of(permutations{:});
            for i = 1:length(self.relators)
                assert(self.niceGroup.isIdentity(self.niceMonomorphismImage(self.relators{i})));
            end
        end

        function res = eqv(self, x, y)
            xyI = self.compose(x, self.inverse(y));
            res = self.niceGroup.isIdentity(self.computeImage(self.niceGroup, xyI));
        end

        function self = FiniteFPGroup(names, relatorLetters, id)
            self@replab.FPGroup(names, relatorLetters, id);
            self.parent = self;
        end

        function p = niceMonomorphismImage(self, w)
            p = w.computeImage(self.niceGroup, self.niceGroup.generators);
        end

        function res = mrdivide(self, rhs)
            if isa(rhs, 'replab.NiceFiniteGroup')
                res = mrdivide@replab.NiceFiniteGroup(self, rhs);
            else
                res = mrdivide@replab.FPGroup(self, rhs);
            end
        end

        function m = morphismByImages(self, target, generatorImages)
            m = morphismByImages@replab.FPGroup(self, target, generatorImages);
        end

    end

end
