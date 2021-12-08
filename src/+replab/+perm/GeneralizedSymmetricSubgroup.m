classdef GeneralizedSymmetricSubgroup < replab.gen.FiniteGroup
% Describes a subgroup of the wreath product of the symmetric group acting on copies of the cyclic group

    methods

        function self = GeneralizedSymmetricSubgroup(n, m, generators, varargin)
        % Constructs a subgroup of the generalized symmetric group
        %
        % Additional keyword arguments (``type``, ``nice`` and ``niceIsomorphism``) are used internally
        % but are not part of the public API.
        %
        % Args:
        %   n (integer): Number of copies of the cyclic group
        %   m (integer): Order of each cyclic group
        %   generators (cell(1,\*) of group elements): Group generators
        %
        % Keyword Args:
        %   generatorNames (cell(1,\*) of charstring): Names of the generators
        %   order (vpi, optional): Order of the group
        %   relators (cell(1,\*) of charstring): Relators given either in word or letter format
            args = struct('type', []);
            [args, rest] = replab.util.populateStruct(args, varargin);
            if isempty(args.type)
                type = replab.perm.GeneralizedSymmetricGroupType(n, m);
            else
                type = args.type;
            end
            self@replab.gen.FiniteGroup(type, generators, rest{:});
            assert(~isempty(self.identity));
        end

        function res = m(self)
            res = self.type.m;
        end

        function res = n(self)
            res = self.type.n;
        end

    end

    methods % Matrix representation of group elements

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

end
