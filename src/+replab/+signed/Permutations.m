classdef Permutations < replab.signed.PermutationGroup
% Describes the signed permutation group over
% {-n...-1, 1...n} where n = domainSize
    
    methods
        
        function self = Permutations(domainSize)
            self.identity = 1:domainSize;
            self.domainSize = domainSize;
            self.parent = self;
            switch self.domainSize
              case 0
                self.generators = cell(1, 0);
              case 1
                self.generators = {[-1]};
              case 2
                self.generators = {[2 1] [-1 2]};
              otherwise
                shift = [2:self.domainSize 1];
                trans = [2 1 3:self.domainSize];
                flip = [-1 2:self.domainSize];
                self.generators = {shift trans flip};
            end
        end
        
        %% Str methods
                
        function s = headerStr(self)
            s = sprintf('Signed permutations acting on {-%d..-1 1..%d}', self.domainSize, self.domainSize);
        end
        
        %% Domain methods
        
        function s = sample(self)
            n = self.domainSize;
            s = randperm(n) .* (randi([0 1], 1, n)*2-1);
        end
        
        %% Finite group methods
                
        function b = contains(self, g)
            assert(length(g) == self.domainSize, 'Signed permutation in wrong domain');
            assert(all(g ~= 0) && all(abs(g) <= self.domainSize), ...
                   'Signed permutation has out of range coefficients');
            b = true;
        end
        
        
        function G = subgroup(self, generators, order)
            if nargin < 3
                order = [];
            end
            G = replab.signed.PermutationSubgroup(self, generators, order);
        end
        
    end
    
    methods (Access = protected)
        
        function o = computeOrder(self)
            o = factorial(vpi(self.domainSize)) * vpi(2)^self.domainSize;
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
            els = [-n:-1 1:n];
            for i = 1:n
                ind0 = ind0 * 2*(n - i + 1);
                ind0 = ind0 + (find(els == g(i)) - 1);
                els = setdiff(els, [g(i) -g(i)]);
            end
            ind = ind0 + 1;
        end
        
        function g = enumeratorAt(self, ind)
            n = self.domainSize;
            ind0 = ind - 1; % make it 0-based
            inds = zeros(1, n);
            for i = 1:n
                r = mod(ind0, 2*i);
                ind0 = (ind0 - r)/(2*i);
                inds(i) = double(r + 1);
            end
            inds = fliplr(inds);
            els = [-n:-1 1:n];
            g = zeros(1, n);
            for i = 1:n
                e = els(inds(i));
                g(i) = e;
                els = setdiff(els, [e -e]);
            end
        end

    end
        
    methods (Static)
        
        function perm = toPermutation(signedPerm)
        % Returns the permutation corresponding to the given signed permutation
        % where the permutation acts on the list [-d,..,-1,1,..,d]
            n = length(signedPerm);
            perm = zeros(1, 2*n);
            for i = 1:length(signedPerm)
                im = n - i + 1; % position of -i in the list
                ip = n + i; % position of i in the list
                j = signedPerm(i);
                jm = n - abs(j) + 1;
                jp = n + abs(j);
                if j > 0
                    perm(im) = jm;
                    perm(ip) = jp;
                else
                    perm(im) = jp;
                    perm(ip) = jm;
                end
            end
        end
        
        function signedPerm = fromPermutation(perm)
        % Returns the signed permutation corresponding to the given permutation encoding
        % see SignedPermutations.toPermutation
            n2 = length(perm);
            if mod(n2, 2) ~= 0
                error('Not an image of a signed permutation');
            end
            n = n2/2; % domain size of the signed permutation
            perm(perm <= n) = -(n - perm(perm <= n) + 1);
            perm(perm > n) = perm(perm > n) - n;
            mperm = -fliplr(perm(1:n));
            pperm = perm(n+1:n2);
            assert(isequal(mperm, pperm), 'Not an image of a signed permutation');
            signedPerm = pperm;
        end

        function mat = toMatrix(signedPerm)
        % Returns the signed permutation matrix corresponding to the given
        % signed permutation such that matrix multiplication is
        % compatible with composition of permutations, i.e. 
        % S.toMatrix(S.compose(x, y)) = 
        % S.toMatrix(x) * S.toMatrix(y)
        % where S = SignedPermutations(domainSize)
            n = length(signedPerm);
            mat = sparse(abs(signedPerm), 1:n, sign(signedPerm), n, n);
            if ~replab.Parameters.useSparse
                mat = full(mat);
            end
        end
        
        function b = isSignedPermutationMatrix(mat)
        % Returns true when "mat" is a signed permutation matrix, i.e. a monomial matrix
        % with nonzero entries equal to +1 or -1
            if isequal(size(mat), [0 0])
                b = true;
                return
            end
            n = size(mat, 1);
            [I J S] = find(mat);
            I = I';
            J = J';
            S = S';
            sI = sort(I);
            [sJ IJ] = sort(J);
            b = isequal(sI, 1:n) && isequal(sJ, 1:n) && isequal(abs(S), ones(1, n));
        end
        
        function signedPerm = fromMatrix(mat)
        % Returns the signed permutation corresponding to the given matrix representation
        % or throws an error
            if isequal(size(mat), [0 0])
                signedPerm = zeros(1, 0);
                return
            end
            signedPerm = [];
            n = size(mat, 1);
            [I J S] = find(mat);
            if length(I) ~= n
                error('Not a monomial matrix');
            end
            I = I';
            J = J';
            S = S';
            sI = sort(I);
            [sJ IJ] = sort(J);
            if ~isequal(sI, 1:n) || ~isequal(sJ, 1:n)
                error('Not a monomial matrix');
            end
            if ~isequal(abs(S), ones(1, n))
                error('Monomial matrix with entries other than +1, -1');
            end
            signedPerm = I.*S;
            signedPerm = signedPerm(IJ);
        end
    
    end
    
    methods (Static)
        
        function Q = quaternionGroup(self)
        % Returns a signed representation of the quaternion group
            SS4 = replab.signed.Permutations(4);
            g1 = [-1 -2 -3 -4];
            gi = [2 -1 4 -3];
            gj = [3 -4 -1 2];
            Q = SS4.subgroup({g1 gi gj});
        end
        
    end
    
end
