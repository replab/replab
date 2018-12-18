classdef Permutations < replab.cat.Group
   
    properties (SetAccess = protected)
        canEqv;
        canHash;
        canSample;
        domainSize;
    end
    
    methods
        
        function self = Permutations(domainSize)
            self.canEqv = true;
            self.canHash = true;
            self.canSample = true;
            self.domainSize = domainSize;
            self.identity = 1:domainSize;
        end
        
        function p = fromCycles(varargin)
        % Constructs a permutation from a product of cycles, each
        % cycle being a row vector, and the sequence cycles being
        % given as variable arguments
            n = self.domainSize;
            for i = length(varargin):-1:1
                cycle = varargin{i};
                % cycle 2 3 1 means that 2 -> 3, 3 -> 1, 1 -> 2
                cycleImage = [cycle(2:end) cycle(1)];
                newEl = 1:n;
                newEl(cycle) = cycleImage;
                p = replab.Perm.compose(newEl, p);
            end
        end
        
        function s = str(self)
            s = sprintf('Permutations acting on %d elements', self.domainSize);
        end
        
        function b = eqv(self, x, y)
            b = isequal(x, y);
        end
        
        function h = hash(self, x)
            h = replab.cat.Domain.hashIntegers(x);
        end
        
        function s = sample(self)
            s = randperm(self.domainSize);
        end
        
        function z = compose(self, x, y)
            z = x(y);
        end
        
        function y = inverse(self, x)
            y = zeros(1, self.domainSize);
            y(x) = 1:self.domainSize;
        end
        
        function A = naturalAction(self)
            n = self.domainSize;
            desc = sprintf('Permutations of size %d acting on integers [1...%d]', n, n);
            P = replab.cat.Domain.integerRange(1, n);
            A = replab.cat.BSGSActionFun(desc, self, P, ...
                                         @(g, p) g(p), ...
                                         @(g) replab.Permutations.findMovedPoint_(g));
        end

        function A = vectorAction(self, field)
            n = self.domainSize;
            desc = sprintf('Permutations acting on %d dimensional vectors in %s', n, field);
            P = replab.Vectors(n, field);
            A = replab.cat.ActionFun(desc, self, P, ...
                                     @(g, p) replab.Permutations.vectorImage_(g, p));
        end
        
        function A = selfAdjointMatrixAction(self, field)
            n = self.domainSize;
            desc = sprintf('Permutations acting on %d x %d self-adjoint matrices in %s', n, field);
            P = replab.SelfAdjointMatrices(n, field);
            A = replab.cat.ActionFun(desc, self, P, ...
                                     @(g, p) replab.Permutations.selfAdjointMatrixImage_(g, p));
        end

        function R = naturalRepresentation(self)
            group = replab.PermutationGroup.symmetric(self.domainSize);
            fun = @(x) replab.Permutations.toMatrix_(x);
            target = replab.PermutationMatrices(self.domainSize);
            R = replab.RepFun(group, fun, self, target);
        end
        
    end
    
    methods (Static, Access = protected)
        
        function p = findMovedPoint_(perm)
            for i = 1:length(perm)
                if i ~= perm(i)
                    p = i;
                    return
                end
            end
            p = [];
        end

        function vec = vectorImage_(perm, vec)
        % Permutation of a column vector
        % equivalent to Perm.matrix(perm) * vec
            assert(length(perm) == length(vec));
            vec(perm) = vec;
        end
        
        function M = selfAdjointMatrixImage_(perm, M)
        % Returns the image under the action of perm on the columns and rows of M
        % i.e. Perm.matrix(perm)*M*Perm.matrix(perm)'
            n = length(perm);
            assert(size(M, 1) == n);
            assert(size(M, 2) == n);
            M(perm, perm) = M;
        end
        
        function mat = toMatrix_(perm)
        % Returns the permutation matrix corresponding to the given permutation
        % such that matrix multiplication is compatible with composition of
        % permutations, i.e. 
        % Perm.matrix(Perm.compose(x, y)) = Perm.matrix(x) * Perm.matrix(y)
            n = length(perm);
            mat = sparse(perm, 1:n, ones(1, n), n, n);
        end

    end

end
