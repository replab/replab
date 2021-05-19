classdef GeneralizedSymmetricGroup < replab.perm.GeneralizedSymmetricSubgroup

    methods

        function self = GeneralizedSymmetricGroup(n, m)
        % Constructs the generalized symmetric group
        %
        % Args:
        %   n (integer): Number of copies of the cyclic group
        %   m (integer): Order of each cyclic group
            o = replab.util.multiplyIntegers(1:n)*replab.util.multiplyIntegers(ones(1,n)*m);
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
            self@replab.perm.GeneralizedSymmetricSubgroup(n, m, generators, 'order', o, 'type', 'self');
        end

        function m = naturalMorphism(self, larger)
        % Returns the natural morphism to a GeneralizedSymmetricGroup
        %
        % The given group needs to be larger/compatible in the sense that ``self.n == larger.n`` and
        % ``self.m`` divides ``larger.m``.
        %
        % Args:
        %   gsg (`+replab.+perm.GeneralizedSymmetricGroup`): Generalized symmetric group in the sense above
        %
        % Returns:
        %   `+replab.FiniteMorphism`: An injective morphism
            assert(self.n == larger.n, 'Groups must have the same size');
            assert(mod(larger.m, self.m) == 0, 'Cyclic orders must be compatible');
            factor = larger.m/self.m;
            m = self.morphismByFunction(larger, @(g) [g(1,:); g(2,:)*factor]);
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

        function [G g] = fromMatrix(M)
        % Recognizes a generalized permutation matrix
        %
        % Returns
        % -------
        %   G:
        %     `.GeneralizedSymmetricGroup` or ``[]``: Group the recognized element is part of, or ``[]`` if unsuccessful
        %   g:
        %     group element or ``[]``: Group element if successful otherwise ``[]``
            G = [];
            g = [];
            [I J V] = find(M);
            n = length(I);
            sI = sort(I);
            [sJ IJ] = sort(J);
            if ~all(sI' == 1:n) || ~all(sJ' == 1:n)
                return
            end
            if any(V.*conj(V) ~= 1)
                return
            end
            angles = angle(double(V));
            angles(angles<0) = angles(angles<0) + 2*pi;
            [N D] = rat(angles/2/pi);
            o = 1;
            for i = 1:n
                if V(i) == 1i
                    o = lcm(o, 4);
                    N(i) = 1;
                    D(i) = 4;
                elseif V(i) == -1i
                    o = lcm(o, 4);
                    N(i) = 3;
                    D(i) = 4;
                else
                    o = lcm(o, D(i));
                    if replab.cyclotomic.E(D(i))^N(i) ~= V(i)
                        return
                    end
                end
            end
            N = N.*o./D;
            G = replab.perm.GeneralizedSymmetricGroup(n, o);
            g = [I(IJ)'; N(IJ)'];
        end

    end

end
