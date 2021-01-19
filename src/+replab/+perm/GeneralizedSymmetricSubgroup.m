classdef GeneralizedSymmetricSubgroup < replab.NiceFiniteGroup
% Describes a subgroup of the wreath product of the symmetric group acting on copies of the cyclic group

    properties
        n % (integer): Number of copies of the cyclic group
        m % (integer): Order of each cyclic group
    end

    methods

        function self = GeneralizedSymmetricSubgroup(n, m, generators, order, type)
        % Constructs a subgroup of the generalized symmetric group
        %
        % Args:
        %   n (integer): Number of copies of the cyclic group
        %   m (integer): Order of each cyclic group
        %   generators (cell(1,\*) of group elements): Group generators
        %   order (vpi, optional): Order of the group
        %   type (`+replab.SignedPermutationGroup`, optional): Type of this group if known, or ``'self'`` if this group is its own type
            self.n = n;
            self.m = m;
            self.identity = [1:n; zeros(1, n)];
            self.generators = generators;
            if nargin > 2 && ~isempty(order)
                self.cache('order', order, '==');
            end
            if nargin < 4
                type = [];
            end
            if isempty(type)
                self.type = replab.perm.GeneralizedSymmetricGroup(n, m);
            elseif isequal(type, 'self')
                self.type = self;
            else
                self.type = type;
            end

        end

    end

    methods

        function rho = naturalRep(self, field)
        % Returns the natural representation of this group using generalized permutation matrices
        %
        % Args:
        %   field ({'R', 'C'}): Whether the representation is real or complex
        %
        % Returns:
        %   `+replab.Rep`: Representation
            rho = replab.rep.GeneralizedPermutationNaturalRep(self, field);
        end

        function M = toCyclotomicMatrix(self, x)
        % Returns the generalized permutation matrix corresponding to the given element
        %
        % Args:
        %   x (group element): Element to compute the matrix representation of
        %
        % Returns:
        %   `+replab.cyclotomic`: Cyclotomic matrix
            M = replab.cyclotomic.zeros(self.n, self.n);
            for i = 1:self.n
                M(x(1,i), i) = replab.cyclotomic.E(self.m)^x(2,i);
            end
        end

        function M = toIntvalMatrix(self, x)
        % Returns the generalized permutation matrix corresponding to the given element
        %
        % Args:
        %   x (group element): Element to compute the matrix representation of
        %
        % Returns:
        %   `+replab.cyclotomic`: Cyclotomic matrix
            M = replab.cyclotomic.zeros(self.n, self.n);
            V = exp(2i*intval('pi')*intval(x(2,:))/intval(self.m));
            V(x(2,:) == 0) = intval(1);
            V(2*x(2,:) == self.m) = intval(-1);
            V(4*x(2,:) == self.m) = intval(1i);
            V(4*x(2,:) == 3*self.m) = intval(-1i);
            M = sparse(x(1,:), 1:self.n, V, self.n, self.n);
        end

        function M = toSparseMatrix(self, x)
        % Returns the generalized permutation matrix corresponding to the given element
        %
        % Args:
        %   x (group element): Element to compute the matrix representation of
        %
        % Returns:
        %   double(\*,\*): Double sparse matrix
            V = exp(2i*pi*x(2,:)/self.m);
            V(x(2,:) == 0) = 1;
            V(2*x(2,:) == self.m) = -1;
            V(4*x(2,:) == self.m) = 1i;
            V(4*x(2,:) == 3*self.m) = -1i;
            M = sparse(x(1,:), 1:self.n, V, self.n, self.n);
        end

        function M = toMatrix(self, x)
        % Returns the generalized permutation matrix corresponding to the given element
        %
        % Args:
        %   x (group element): Element to compute the matrix representation of
        %
        % Returns:
        %   double(\*,\*): Double matrix
            M = full(self.toSparseMatrix(x));
        end

    end

    methods % Implementations

        % Domain

        function b = eqv(self, x, y)
            b = all(all(x == y));
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

        % NiceFiniteGroup methods

        function res = hasSameTypeAs(self, rhs)
            res = isa(rhs, 'replab.perm.GeneralizedSymmetricSubgroup') && (self.type.n == rhs.type.n) && (self.type.m == rhs.type.m);
        end

        function p = nicePreimage(self, p1)
            p = zeros(2, self.n);
            for i = 1:self.n
                v = p1(self.m*(i-1)+1) - 1;
                p(2,i) = mod(v, self.m);
                p(1,i) = (v - p(2,i))/self.m + 1;
            end
        end

        function p1 = niceImage(self, p)
            p1 = zeros(self.m, self.n);
            for i = 1:self.n
                p1(:,i) = mod((0:self.m-1)+p(2,i), self.m) + self.m*(p(1,i)-1) + 1;
            end
            p1 = p1(:)';
        end

        function grp = niceSubgroup(self, generators, order, niceGroup)
        % Constructs a permutation subgroup from its generators
        %
        % Args:
        %   generators (cell(1,\*) of generalized permutations): List of generators given as a permutations in a row cell array
        %   order (vpi, optional): Argument specifying the group order, if given can speed up computations
        %   niceGroup (`+replab.PermutationGroup`, optional): Image of this subgroup under the nice morphism
        %
        % Returns:
        %   `+replab.+perm.GeneralizedSymmetricSubgroup`: The constructed signed permutation subgroup
            if nargin < 4
                niceGroup = [];
            end
            if nargin < 3
                order = [];
            end
            grp = replab.perm.GeneralizedSymmetricSubgroup(self.n, self.m, generators, order, self.type);
            if ~isempty(niceGroup)
                grp.cache('niceGroup', niceGroup, '==');
            end
        end

    end

end
