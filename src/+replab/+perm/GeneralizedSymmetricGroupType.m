classdef GeneralizedSymmetricGroupType < replab.gen.StaticFiniteGroupType
% Wreath product of the symmetric group acting on copies of the cyclic group

    properties (SetAccess = protected)
        n % (integer): Number of copies of the cyclic group
        m % (integer): Order of each cyclic group
    end

    methods

        function self = GeneralizedSymmetricGroupType(n, m)
            self.n = n;
            self.m = m;
            orderFun = @() replab.util.factorial(n)*replab.util.multiplyIntegers(ones(1,n)*m);
            if n < 2
                sourceGenerators = cell(1, 0);
            elseif n == 2
                sourceGenerators = {[2 1; 0 0]};
            else
                sourceGenerators = {[2:n 1; zeros(1, n)] [2 1 3:n; zeros(1, n)]};
            end
            if m > 1
                sourceGenerators{1,end+1} = [1:n; 1 zeros(1, n-1)];
            end
            sourceArgs = {'order', orderFun};
            targetType = replab.perm.PermutationGroupType.make(n*m);
            self.identity = [1:n; zeros(1, n)];
            self.finishConstruction(sourceGenerators, sourceArgs, targetType);
        end

        function mu = naturalMorphismTo(self, targetType)
        % Returns the natural morphism between generalized symmetric groups
        %
        % Both the source (self) and the target group have the same number of copies of the cyclic
        % subgroup (``sourceType.n == targetType.n``).
        % The target cyclic group embeds the source cyclic group in the sense
        % that ``sourceType.m`` divides ``targetType.m``.
        %
        % The method returns a morphism between ``sourceType.parentGroup`` and
        % ``targetType.parentGroup``.
        %
        % Args:
        %   targetType (`+replab.+perm.GeneralizedSymmetricGroupType`): Target type
        %
        % Returns:
        %   `+replab.FiniteMorphism`: An injective morphism
            sourceType = self;
            assert(sourceType.n == targetType.n, 'Number of copies of the cyclic group must be the same');
            assert(mod(targetType.m, sourceType.m) == 0, 'Cyclic orders must be compatible');
            factor = targetType.m/sourceType.m;
            mu = sourceType.parentGroup.morphismByFunction(targetType.parentGroup, @(g) [g(1,:); g(2,:)*factor]);
        end


    end

    methods % Implementations

        % Domain

        function b = eqv(self, x, y)
            b = all(all(x == y));
        end

        function s = sample(self)
            s = [randperm(self.n); randi([0 self.m-1], 1, self.n)];
        end

        % Monoid

        function z = compose(self, x, y)
            z = zeros(2, self.n);
            z(1,:) = x(1,y(1,:));
            z(2,:) = mod(x(2,y(1,:)) + y(2,:), self.m);
        end

        % Group

        function y = inverse(self, x)
            n = self.n;
            m = self.m;
            y = zeros(2, n);
            y(1,x(1,:)) = 1:n;
            y(2,x(1,:)) = mod(self.m - x(2,:), self.m);
        end

        % FiniteGroupType

        function l = isSameTypeAs(self, otherType)
            l = isa(otherType, 'replab.perm.GeneralizedSymmetricGroupType') && self.n == otherType.n && self.m == otherType.m;
        end

        % StaticFiniteGroupType

        function p = preimageElement(self, p1)
            p = zeros(2, self.n);
            for i = 1:self.n
                v = p1(self.m*(i-1)+1) - 1;
                p(2,i) = mod(v, self.m);
                p(1,i) = (v - p(2,i))/self.m + 1;
            end
        end

        function p1 = imageElement(self, p)
            p1 = zeros(self.m, self.n);
            for i = 1:self.n
                p1(:,i) = mod((0:self.m-1)+p(2,i), self.m) + self.m*(p(1,i)-1) + 1;
            end
            p1 = p1(:)';
        end

    end

    methods (Access = protected)

        function G = groupFromNiceImage_(self, generators, nice, niceIsomorphism)

            G = replab.perm.GeneralizedSymmetricSubgroup(self.n, self.m, generators, 'type', self, 'nice', nice, 'niceIsomorphism', niceIsomorphism);
        end

    end

    methods (Static)

        function mu = isomorphismPermutationGroup(G)
        % Returns the isomorphism between a permutation group and a generalized symmetric group
            targetType = replab.perm.GeneralizedSymmetricGroupType(G.domainSize, 1);
            mu = replab.Morphism.lambda(G, targetType.parentGroup, @(g) [g; zeros(1, G.domainSize)]);
        end


        function mu = isomorphismSignedPermutationGroup(G)
            targetType = replab.perm.GeneralizedSymmetricGroupType(G.domainSize, 2);
            mu = replab.Morphism.lambda(G, targetType.parentGroup, @(g) [abs(g); (1-sign(g))/2]);
        end

        function [G, g] = fromMatrix(M)
        % Recognizes a generalized permutation matrix
        %
        % Returns
        % -------
        %   G:
        %     `.GeneralizedSymmetricGroupType` or ``[]``: Group the recognized element is part of, or ``[]`` if unsuccessful
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
                if V(i) == 1
                    N(i) = 1;
                    D(i) = 1;
                elseif V(i) == -1
                    o = lcm(o, 2);
                    N(i) = 1;
                    D(i) = 2;
                elseif V(i) == 1i
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
            G = replab.perm.GeneralizedSymmetricGroupType(n, o);
            g = [I(IJ)'; N(IJ)'];
        end

    end

end
