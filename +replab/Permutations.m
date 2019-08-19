classdef Permutations < replab.PermutationGroup & replab.FiniteGroup
% Describes permutations over n = "domainSize" elements, i.e.
% the symmetric group Sn
    
    methods % Implementations of abstract methods
        
        function self = Permutations(domainSize)
            self@replab.PermutationGroup(domainSize);
            if self.domainSize < 2
                self.generators = cell(1, 0);
            elseif self.domainSize == 2
                self.generators = {[2 1]};
            else
                self.generators = {[2:domainSize 1] [2 1 3:domainSize]};
            end
        end
        
        % Str
                
        function s = headerStr(self)
            s = sprintf('Permutations acting on %d elements', self.domainSize);
        end

        % Domain
        
        function s = sample(self)
            s = randperm(self.domainSize);
        end
        
        % FinitelyGeneratedGroup
        
        function w = factorization(self, x)
        % Factorizes a permutation using bubble sort
            if self.isIdentity(x)
                w = replab.Word.identity;
                return
            elseif self.domainSize == 2
                % not identity
                w = replab.Word.generator(1);
                return
            end
            n = length(x);
            w = replab.Word.identity;
            moved = true;
            while moved
                moved = false;
                for i = 1:n-1
                    if x(i) > x(i+1)
                        t = x(i+1);
                        x(i+1) = x(i);
                        x(i) = t;
                        moved = true;
                        if i == 1
                            shift = replab.Word.identity;
                        else
                            shift = replab.Word.fromIndicesAndExponents(1, i - 1);
                        end
                        w = shift * replab.Word.generator(2) * inv(shift) * w;
                    end
                end
            end
        end
        
        % FiniteGroup
        
        function b = contains(self, g)
            b = (length(g) == self.domainSize) && all(g > 0);
        end
        
        function b = knownOrder(self)
            b = true;
        end
        
        function o = order(self)
            o = factorial(vpi(self.domainSize));
        end
        
        function E = elements(self)
            E = replab.Enumerator.lambda(self.order, ...
                                         @(ind) self.enumeratorAt(ind), ...
                                         @(el) self.enumeratorFind(el));
        end
        
        function d = decomposition(self)
            G = self.subgroup(self.generators, self.order);
            d = G.decomposition;
        end
        
    end
    
    methods (Access = protected)
        
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
        %   i: Shift so that $j$ is sent to $j+i$ (wrapping around).
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

        function grp = subgroup(self, generators, orderOpt)
        % Constructs a permutation subgroup from its generators
        %
        % Args:
        %   generators: List of generators given as a permutations in a row cell array
        %   orderOpt: Optional argument specifying the group order, will speed up computations
        %
        % Returns:
        %   +replab.PermutationSubgroup: The constructed permutation subgroup.
            if nargin < 3
                orderOpt = [];
            end
            grp = replab.perm.PermutationBSGSGroup(self, generators, orderOpt);
        end
        
        function grp = trivialSubgroup(self)
            grp = self.subgroup({}, vpi(1));
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
        % (can be seen as the multiplication of quaternions +/- 1,i,j,k = 8 elements)
            S8 = replab.Permutations(8);
            g1 = S8.fromCycles([1 2 4 7], [3 6 8 5]);
            g2 = S8.fromCycles([1 3 4 8], [2 5 7 6]);
            Q = replab.Permutations(8).subgroup({g1 g2});
        end
        
        function mat = toMatrix(perm)
        % Returns the permutation matrix corresponding to the given permutation
        %
        % The returned matrix is such that matrix multiplication is compatible with composition of
        % permutations, i.e. for `P = replab.Permutations(domainSize)` we have
        % `P.toMatrix(P.compose(x, y)) = P.toMatrix(x) * P.toMatrix(y)`
        %
        % Args:
        %   perm: Permutation
        %
        % Returns:
        %   The permutation matrix corresponding to ``perm``.
            n = length(perm);
            mat = sparse(perm, 1:n, ones(1, n), n, n);
        end
        
        function perm = fromMatrix(mat)
        % Returns the signed permutation corresponding to the given matrix representation
        %
        % See :method:`+replab.Permutations.toMatrix`
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
