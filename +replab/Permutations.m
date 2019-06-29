classdef Permutations < replab.PermutationGroup & replab.FiniteGroup
% Describes permutations over n = "domainSize" elements, i.e.
% the symmetric group Sn
    
    methods % Implementations of abstract methods
        
        function self = Permutations(domainSize)
            self = self@replab.PermutationGroup(domainSize);
            if self.domainSize < 2
                self.generators = cell(1, 0);
            elseif self.domainSize == 2
                self.generators = {[2 1]};
            else
                self.generators = {[2:domainSize 1] [2 1 3:domainSize]};
            end
        end
        
        % Str
                
        function s = shortStr(self, maxColumns)
            s = sprintf('Permutations acting on %d elements', self.domainSize);
        end

        function lines = longStr(self, maxRows, maxColumns)
            lines = replab.str.longStr(self, maxRows, maxColumns);
            lines{1} = self.shortStr(maxColumns);
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
            E = replab.EnumeratorFun(self.order, ...
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
        
        function p = fromCycles(self, varargin)
        % Constructs a permutation from a product of cycles, each
        % cycle being a row vector, and the sequence cycles being
        % given as variable arguments
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
            if nargin < 3
                orderOpt = [];
            end
            grp = replab.perm.PermutationBSGSGroup(self, generators, orderOpt);
        end
        
        function grp = trivialGroup(self)
            grp = self.subgroup({}, vpi(1));
        end
        
        function grp = cyclicGroup(self)
            n = self.domainSize;
            if n == 1
                grp = self.trivialGroup;
            else
                grp = self.subgroup({[2:n 1]}, vpi(n));
            end
        end
        
        function grp = alternatingGroup(self)
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
        
    end
    
    methods (Static)
        
        function mat = toMatrix(perm)
        % Returns the permutation matrix corresponding to the given permutation
        % such that matrix multiplication is compatible with composition of
        % permutations, i.e. for P = replab.Permutations(domainSize)
        % P.toMatrix(P.compose(x, y)) = P.toMatrix(x) * P.toMatrix(y)
            n = length(perm);
            mat = sparse(perm, 1:n, ones(1, n), n, n);
        end
        
        function perm = fromMatrix(mat)
        % Returns the signed permutation corresponding to the given matrix representation
        % or throws an error
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
