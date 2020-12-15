classdef WreathProductGroup_finite < replab.WreathProductGroup & replab.prods.SemidirectProductGroup_finite
% Wreath product of a permutation group acting on a finite group

    methods

        function self = WreathProductGroup_finite(H, A)
            assert(isa(H, 'replab.PermutationGroup'));
            n = H.domainSize;
            base = A.directPower(n);
            phi = replab.perm.PermutationCellAction(H, base);
            self@replab.prods.SemidirectProductGroup_finite(phi);
            self.n = n;
            self.A = A;
        end

    end

    methods % Implementation

        % Domain

        function g = sample(self)
            g = sample@replab.prods.SemidirectProductGroup_finite(self); % force method selection
        end

        function b = eqv(self, x, y)
            b = eqv@replab.prods.SemidirectProductGroup_finite(self, x, y);
        end

        % Monoid

        function z = compose(self, x, y)
            z = compose@replab.prods.SemidirectProductGroup_finite(self, x, y);
        end

        % Group

        function xInv = inverse(self, x)
            xInv = inverse@replab.prods.SemidirectProductGroup_finite(self, x);
        end

        % NiceFiniteGroup

        function p = niceImage(self, w)
            p = self.imprimitivePermutation(w, @(a) self.A.niceMorphism.imageElement(a));
        end

        function res = hasSameTypeAs(self, rhs)
            res = isa(rhs, 'replab.WreathProductGroup') && self.n == rhs.n && self.A.hasSameTypeAs(rhs.A);
        end

    end

end
