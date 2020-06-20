classdef Permutations < replab.PermutationGroup
% Describes permutations over n = "domainSize" elements, i.e. the symmetric group Sn
%
% Example:
%   >>> S5 = replab.Permutations(5);
%   >>> S5.order
%      ans =
%      120

    methods % Implementations of abstract methods

        function self = Permutations(domainSize)
        % Constructs the symmetric over a given domain size
        %
        % Args:
        %   domainSize (integer): Domain size, must be > 0
            o = factorial(vpi(domainSize));
            if domainSize < 2
                generators = cell(1, 0);
            elseif domainSize == 2
                generators = {[2 1]};
            else
                generators = {[2:domainSize 1] [2 1 3:domainSize]};
            end
            self = self@replab.PermutationGroup(domainSize, generators, o, []);
        end

        %% Str methods

        function s = headerStr(self)
            s = sprintf('Permutations acting on %d elements', self.domainSize);
        end

        %% Domain methods

        function s = sample(self)
            s = randperm(self.domainSize); % overriden for efficiency
        end

        %% FiniteGroup methods

        function b = contains(self, g)
            assert(length(g) == self.domainSize, 'Permutation in wrong domain');
            assert(all(g > 0), 'Permutation should have positive coefficients');
            b = true;
        end

    end

    methods (Access = protected)

        function o = computeOrder(self)
            o = factorial(vpi(self.domainSize));
        end

        function E = computeElements(self)
            E = replab.IndexedFamily.lambda(self.order, ...
                                            @(ind) self.enumeratorAt(ind), ...
                                            @(el) self.enumeratorFind(el));
        end

        function d = computeDecomposition(self)
            G = self.subgroup(self.generators, self.order);
            d = G.decomposition;
        end

        function ind = enumeratorFind(self, g)
            n = self.domainSize;
            ind0 = vpi(0);
            els = 1:n;
            for i = 1:n
                ind0 = ind0 * (n - i + 1);
                ind0 = ind0 + (find(els == g(i)) - 1);
                els = setdiff(els, g(i));
            end
            ind = ind0 + 1;
        end

        function g = enumeratorAt(self, ind)
            n = self.domainSize;
            ind0 = ind - 1; % make it 0-based
            inds = zeros(1, n);
            for i = 1:n
                r = mod(ind0, i);
                ind0 = (ind0 - r)/i;
                inds(i) = double(r + 1);
            end
            inds = fliplr(inds);
            els = 1:n;
            g = zeros(1, n);
            for i = 1:n
                e = els(inds(i));
                g(i) = e;
                els = setdiff(els, e);
            end
        end

    end

    methods

        function p = transposition(self, i, j)
        % Returns the transposition permuting ``i`` and ``j``.
        %
        % Args:
        %   i: First domain element to be transposed.
        %   j: Second domain element to be transposed.
        %
        % Returns:
        %   The constructed transposition.
            assert(1 <= i);
            assert(i <= self.domainSize);
            assert(1 <= j);
            assert(j <= self.domainSize);
            assert(i ~= j);
            p = 1:self.domainSize;
            p([i j]) = [j i];
        end

        function p = shift(self, i)
        % Returns the cyclic permutation that shifts the domain indices by ``i``.
        %
        % Args:
        %   i: Shift so that ``j`` is sent to ``j + i`` (wrapping around).
        %
        % Returns:
        %   The constructed cyclic shift.
            n = self.domainSize;
            p = mod((0:n-1)+i, n)+1;
        end

        function p = fromCycles(self, varargin)
        % Constructs a permutation from a product of cycles.
        %
        % Each cycle is given as a row vector, and the sequence of cycles is
        % given as variable arguments.
        %
        % Args:
        %  *args: Sequence of cycles as row vectors of indices
        %
        % Returns:
        %  The permutation corresponding to the product of cycles.
            n = self.domainSize;
            p = self.identity;
            for i = length(varargin):-1:1
                cycle = varargin{i};
                % cycle 2 3 1 means that 2 -> 3, 3 -> 1, 1 -> 2
                cycleImage = [cycle(2:end) cycle(1)];
                newEl = 1:n;
                newEl(cycle) = cycleImage;
                p = self.compose(newEl, p);
            end
        end

        function grp = cyclicSubgroup(self)
            n = self.domainSize;
            if n == 1
                grp = self.trivialGroup;
            else
                grp = self.subgroup({[2:n 1]}, vpi(n));
            end
        end

        function grp = alternatingSubgroup(self)
            n = self.domainSize;
            if n <= 2
                grp = self.trivialGroup;
            else
                c3 = [2 3 1 4:n];
                if mod(n, 2) == 0
                    s = [1 3:n 2];
                else
                    s = [2:n 1];
                end
                grp = self.subgroup({c3 s}, self.order/2);
            end
        end

        function grp = symmetricGroup(self)
            grp = self.subgroup(self.generators, self.order);
        end

        function grp = dihedralSubgroup(self)
            n = self.domainSize;
            if n <= 2
                grp = self.symmetricGroup;
            else
                g1 = [2:n 1];
                g2 = fliplr(1:n);
                grp = self.subgroup({g1 g2});
            end
        end

    end

    methods (Static)

        function Q = quaternionGroup(self)
        % Returns a permutation representation of the quaternion group
        %
        % The quaternion group returned can be seen as the multiplication of
        % of the four quaternion generators 1, i, j, k with a sign +/-1, thus
        % can be represented by permutations on 8 elements.
        %
        % Returns:
        %    +replab.PermutationGroup: Quaternion group as a subgroup of S(8)
            S8 = replab.Permutations(8);
            g1 = S8.fromCycles([1 2 4 7], [3 6 8 5]);
            g2 = S8.fromCycles([1 3 4 8], [2 5 7 6]);
            Q = replab.Permutations(8).subgroup({g1 g2});
        end

        function mat = toSparseMatrix(perm)
        % Returns the sparse permutation matrix corresponding to the given permutation
        %
        % The returned matrix is such that matrix multiplication is compatible with composition of
        % permutations, i.e. for ``P = replab.Permutations(domainSize)`` we have
        % ``P.toMatrix(P.compose(x, y)) = P.toMatrix(x) * P.toMatrix(y)``
        %
        % Args:
        %   perm (permutation row vector): Permutation
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
        %   perm (permutation row vector): Permutation
        %
        % Returns:
        %   The permutation matrix corresponding to ``perm``.
            mat = full(replab.Permutations.toSparseMatrix(perm));
        end

        function perm = fromMatrix(mat)
        % Returns the signed permutation corresponding to the given matrix representation
        %
        % See `+replab.Permutations.toMatrix`
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

    end

end
