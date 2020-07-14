classdef OfFiniteGroup < replab.semidirectproduct.OfFiniteGroups & replab.wreathproduct.Common
% Wreath product of a permutation group acting on a finite group

    methods

        function self = OfFiniteGroup(H, A)
            assert(isa(H, 'replab.PermutationGroup'));
            n = H.domainSize;
            base = A.directPower(n);
            phi = replab.perm.PermutationCellAction(H, base);
            self@replab.semidirectproduct.OfFiniteGroups(phi);
            assert(isa(A, self.requiredType));
            self.n = n;
            self.A = A;
        end

        function t = requiredType(self)
            t = 'replab.FiniteGroup';
        end

    end

    methods % Implementation

        % Domain

        function g = sample(self)
            g = sample@replab.semidirectproduct.OfCompactGroups(self); % force method selection
        end

        function b = eqv(self, x, y)
            b = eqv@replab.semidirectproduct.OfCompactGroups(self, x, y);
        end

        % Monoid

        function z = compose(self, x, y)
            z = compose@replab.semidirectproduct.OfCompactGroups(self, x, y);
        end

        % Group

        function xInv = inverse(self, x)
            xInv = inverse@replab.semidirectproduct.OfCompactGroups(self, x);
        end

        % NiceFiniteGroup

        function p = niceImage(self, w)
            p = self.imprimitivePermutation(w, @(a) self.A.niceMorphism.imageElement(a));
        end

        function res = hasSameTypeAs(self, rhs)
            res = isa(rhs, 'replab.wreathproduct.OfFiniteGroup') && self.n == rhs.n && self.A.hasSameTypeAs(rhs.A);
        end

    end

end
