classdef GeneralizedSymmetricGroup < replab.perm.GeneralizedSymmetricSubgroup

    methods

        function self = GeneralizedSymmetricGroup(n, m)
        % Constructs the generalized symmetric group
        %
        % Args:
        %   n (integer): Number of copies of the cyclic group
        %   m (integer): Order of each cyclic group
            o = factorial(vpi(n))*vpi(m)^n;
            if n < 2
                generators = cell(1, 0);
            elseif n == 2
                generators = {[2 1; 0 0]};
            else
                generators = {[2:n 1; zeros(1, n)] [2 1 3:n; zeros(1, n)]};
            end
            if m > 1
                generators{1,end+1} = [1:n; 1 zeros(1, n-1)];
            end
            self = self@replab.perm.GeneralizedSymmetricSubgroup(n, m, generators, o, 'self');
        end

    end

    methods % Implementations

        % Domain

        function s = sample(self)
            s = [randperm(self.n); randi([0 self.m-1], 1, self.n)];
        end

        % FiniteGroup

        function b = contains(self, g)
            assert(size(g, 2) == self.n, 'Wrong domain size');
            b = true;
        end

    end

    methods (Static)

        function f = morphismFromPermutationGroup(G)
            f = replab.Morphism.lambda(G, replab.perm.GeneralizedSymmetricGroup(G.domainSize, 1), @(g) [g; zeros(1, G.domainSize)]);
        end


        function f = morphismFromSignedPermutationGroup(G)
            f = replab.Morphism.lambda(G, replab.perm.GeneralizedSymmetricGroup(G.domainSize, 2), @(g) [abs(g); (1-sign(g))/2]);
        end

    end

end
