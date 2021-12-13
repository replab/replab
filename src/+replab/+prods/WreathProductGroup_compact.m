classdef WreathProductGroup_compact < replab.WreathProductGroup & replab.prods.SemidirectProductGroup_compact
% Wreath product of a permutation group acting on a compact group

    methods

        function self = WreathProductGroup_compact(H, A)
            assert(isa(H, 'replab.PermutationGroup'));
            n = H.domainSize;
            base = A.directPower(n);
            phi = replab.perm.PermutationCellAction(H, base);
            self@replab.prods.SemidirectProductGroup_compact(phi);
            self.n = n;
            self.A = A;
        end

    end

    methods % Implementations

        % Domain

        function b = eqv(self, x, y)
            b = eqv@replab.SemidirectProductGroup(self, x, y);
        end

        function g = sample(self)
            g = sample@replab.SemidirectProductGroup(self);
        end

        % Monoid

        function z = compose(self, x, y)
            z = compose@replab.SemidirectProductGroup(self, x, y);
        end

        % Group

        function z = inverse(self, x)
            z = inverse@replab.SemidirectProductGroup(self, x);
        end

    end

end
