classdef Permutation
% This class regroups functions that are useful when dealing with permutations.
%
% In RepLAB, permutations are represented using row vectors of integers, those integers being
% encoded using the type ``double``; so this class does not have any instances.

    methods (Static)

        function p = fromCycles(n, varargin)
        % Constructs a permutation from a product of cycles.
        %
        % Each cycle is given as a row vector, and the sequence of cycles is given as variable arguments.
        %
        % Args:
        %   n (integer): Domain size
        %   varargin (cell(1,\*) of integer(1,\*)): Sequence of cycles as row vectors of indices
        %
        % Returns:
        %   permutation: The permutation corresponding to the product of cycles.
            p = 1:n;
            for i = length(varargin):-1:1
                cycle = varargin{i};
                % cycle 2 3 1 means that 2 -> 3, 3 -> 1, 1 -> 2
                cycleImage = [cycle(2:end) cycle(1)];
                newEl = 1:n;
                newEl(cycle) = cycleImage;
                p = newEl(p); % compose(newEl, p);
            end
        end

        function p = transposition(n, i, j)
        % Returns the transposition permuting ``i`` and ``j``.
        %
        % Args:
        %   n (integer): Domain size
        %   i (integer): First domain element to be transposed.
        %   j (integer): Second domain element to be transposed.
        %
        % Returns:
        %   permutation: The constructed transposition
            assert(1 <= i && i <= n);
            assert(1 <= j && j <= n);
            assert(i ~= j);
            p = 1:n;
            p([i j]) = [j i];
        end

        function p = shift(n, i)
        % Returns the cyclic permutation that shifts the domain indices by ``i``.
        %
        % Args:
        %   n (integer): Domain size
        %   i (nonnegative integer): Shift so that ``j`` is sent to ``j + i`` (wrapping around).
        %
        % Returns:
        %   permutation: The constructed cyclic shift
            p = mod((0:n-1)+i, n)+1;
        end

        function o = order(perm)
        % Returns the order of the permutation
        %
        % Args:
        %   perm (permutation): Permutation
        %
        % Returns:
        %   integer: The order of ``perm``, i.e. the smallest ``o`` such that ``perm^o == identity``
            if length(perm) < 2
                o = 1;
                return
            end
            orders = unique(replab.Permutation.cycleStructure(perm));
            if isempty(orders)
                o = 1;
                return
            end
            o = orders(1);
            i = 2;
            for i = 2:length(orders)
                assert(log2(o) + log2(orders(i)) < 53, 'Order of element too big to fit in a double');
                o = lcm(o, orders(i));
            end
        end

        function cycles = cycles(perm)
        % Returns the cycles of the given permutation
        %
        % Example:
        %   >>> isequal(replab.Permutation.cycleStructure([2 1 4 3]), [2 2])
        %       1
        %
        % Args:
        %   perm (permutation): Permutation
        %
        % Returns:
        %   integer(1,\*): Nonincreasing vector of integers containing the sizes of cycles (for sizes > 1)
            x = perm;
            n = length(x);
            cycles = {};
            for i = 1:n
                if x(i) == 0 || x(i) == i % skip 1-cycles or visited points
                    x(i) = 0;
                    continue
                end
                j = i;
                cycle = [];
                while x(j) ~= 0
                    j1 = x(j);
                    x(j) = 0;
                    j = j1;
                    cycle = [cycle j];
                end
                cycles{1,end+1} = cycle;
            end
        end

        function c = cycleStructure(perm)
        % Returns the cycle structure of the given permutation
        %
        % Example:
        %   >>> isequal(replab.Permutation.cycleStructure([2 1 4 3]), [2 2])
        %       1
        %
        % Args:
        %   perm (permutation): Permutation
        %
        % Returns:
        %   integer(1,\*): Nonincreasing vector of integers containing the sizes of cycles (for sizes > 1)
            x = perm;
            n = length(x);
            c = [];
            for i = 1:n
                if x(i) == 0 || x(i) == i % skip 1-cycles or visited points
                    x(i) = 0;
                    continue
                end
                cycleSize = 0;
                j = i;
                while x(j) ~= 0
                    pHold = x(j);
                    x(j) = 0;
                    j = pHold;
                    cycleSize = cycleSize + 1;
                end
                c = [c cycleSize];
            end
            c = fliplr(sort(c));
        end

        function s = sign(perm)
        % Returns the sign of a given permutation
        %
        % Example:
        %   >>> replab.Permutation.sign([3 2 1 4])
        %       -1
        % Args:
        %   perm (permutation): Permutation to compute the sign of
        %
        % Returns:
        %   integer: Sign of the permutation
            x = perm;
            n = length(x);
            s = 1; % Permutation sign, 1 if even, -1 if odd
            for i = 1:n
                if x(i) == 0 || x(i) == i % Skip over 1-cycles and numbers that have been cycled through
                    continue
                end
                cycleSize = 0;
                j = i;
                while x(j) ~= 0
                    pHold = x(j);
                    x(j) = 0;
                    j = pHold;
                    cycleSize = cycleSize + 1;
                end
                if mod(cycleSize, 2) == 0 % cycle size is even
                    s = -s; % flip
                end
            end
        end

        function p = sorting(array, greaterFun)
        % Returns the permutation that sorts a cell array using a custom comparison function
        %
        % Args:
        %   array: A (row or column) vector array containing data
        %          of arbitrary type
        %   greaterFun: A function handle that compares elements such that
        %               ``greaterFun(x, y) == true`` when x > y
        %               and false otherwise
        % Returns:
        %   A permutation ``p`` such that ``sorted = array(p)``
            n = length(array);
            p = 1:n;
            inc = round(n/2);
            while inc > 0
                for i = (inc+1):n
                    tmp = p(i);
                    j = i;
                    if isa(array, 'cell')
                        while (j >= inc+1) && greaterFun(array{p(j-inc)}, array{tmp})
                            p(j) = p(j-inc);
                            j = j - inc;
                        end
                    else
                        while (j >= inc+1) && greaterFun(array(p(j-inc)), array(tmp))
                            p(j) = p(j-inc);
                            j = j - inc;
                        end
                    end
                    p(j) = tmp;
                end
                if inc == 2
                    inc = 1;
                else
                    inc = round(inc/2.2);
                end
            end
        end

        function mat = toSparseMatrix(perm)
        % Returns the sparse permutation matrix corresponding to the given permutation
        %
        % The returned matrix is such that matrix multiplication is compatible with composition of
        % permutations, i.e. for ``P = replab.SymmetricGroup(domainSize)`` we have
        % ``P.toMatrix(P.compose(x, y)) = P.toMatrix(x) * P.toMatrix(y)``
        %
        % Args:
        %   perm (permutation): Permutation
        %
        % Returns:
        %   The sparse permutation matrix corresponding to ``perm``.
            n = length(perm);
            mat = sparse(perm, 1:n, ones(1, n), n, n);
        end

        function mat = toMatrix(perm)
        % Returns the permutation matrix corresponding to the given permutation
        %
        % The returned matrix is such that matrix multiplication is compatible with composition of
        % permutations, i.e. for ``P = replab.Permutations(domainSize)`` we have
        % ``P.toMatrix(P.compose(x, y)) = P.toMatrix(x) * P.toMatrix(y)``
        %
        % Args:
        %   perm (permutation): Permutation
        %
        % Returns:
        %   The permutation matrix corresponding to ``perm``.
            mat = full(replab.Permutation.toSparseMatrix(perm));
        end

        function perm = fromMatrix(mat)
        % Returns the permutation corresponding to the given matrix representation
        %
        % See `+replab.Permutation.toMatrix`
        %
        % Args:
        %   mat: A permutation matrix.
        %
        % Returns:
        %   The permutation corresponding to matrix ``mat``.
        %
        % Raises:
        %   Error: if ``mat`` is not a permutation matrix, throws an error
            if isequal(size(mat), [0 0])
                perm = zeros(1, 0);
                return
            end
            perm = [];
            n = size(mat, 1);
            [I J V] = find(mat);
            if length(I) ~= n
                error('Not a monomial matrix');
            end
            I = I';
            J = J';
            V = V';
            if ~isequal(V, ones(1, n))
                error('Not a permutation matrix');
            end
            sI = sort(I);
            [sJ IJ] = sort(J);
            if ~isequal(sI, 1:n) || ~isequal(sJ, 1:n)
                error('Not a monomial matrix');
            end
            perm = I(IJ);
        end

    end

end
