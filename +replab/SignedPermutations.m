classdef SignedPermutations < replab.cat.Group
   
    properties
        n;
    end
    
    methods
        
        function self = SignedPermutations(n)
            self.n = n;
            self.parentOption = [];
            self.identity = 1:n;
        end
        
        function s = str(self)
            s = sprintf('Signed permutations acting on [-%d..1, 1..%d]', self.n, self.n);
        end
        
        function b = eqv(self, x, y)
            b = isequal(x, y);
        end
        
        function h = hash(self, x)
            h = replab.cat.Domain.hashIntegers(x);
        end
        
        function s = sample(self)
            s = randperm(self.n) .* (randi([0 1], 1, self.n)*2-1);
        end
        
        function z = compose(self, x, y)
            z = x(abs(y)).*sign(y);
        end
        
        function y = inverse(self, x)
            y = zeros(1, self.n);
            xAbs = abs(x);
            y(xAbs) = 1:self.n;
            invFlip = xAbs(x < 0);
            y(invFlip) = -y(invFlip);
        end
        
        function A = naturalAction(self)
            n = self.n;
            desc = sprintf('Signed permutations acting on integers [-%d..-1, 1..%d]', n, n);
            P = replab.cat.Domain.integerRange(1, n);
            A = replab.cat.BSGSActionFun(desc, self, P, ...
                                         @(g, p) g(abs(p))*sign(p), ...
                                         @(g) replab.SignedPermutations.findMovedPoint_(g));
        end
        
        function A = vectorAction(self, field)
            n = self.n;
            desc = sprintf('Signed permutations acting on %d dimensional vectors in %s', n, field);
            P = replab.Vectors(n, field);
            A = replab.cat.ActionFun(desc, self, P, ...
                                     @(g, p) replab.SignedPermutations.vectorImage_(g, p));
        end
        
        function A = selfAdjointMatrixAction(self, field)
            n = self.n;
            desc = sprintf('Signed permutations acting on %d x %d self-adjoint matrices in %s', n, field);
            P = replab.SelfAdjointMatrices(n, field);
            A = replab.cat.ActionFun(desc, self, P, ...
                                     @(g, p) replab.SignedPermutations.selfAdjointMatrixImage_(g, p));
        end
        
        function phi = permutationIsomorphism(self)
            fun = @(x) replab.SignedPermutations.toPermutation_(x);
            target = replab.Permutations(2*self.n);
            phi = replab.cat.GroupMorphismFun(fun, self, target);
        end
        
        function R = naturalRepresentation(self)
            fun = @(x) replab.SignedPermutations.toMatrix_(x);
            target = replab.SignedPermutationMatrices(self.n);
            R = replab.RepFun(fun, self, target);
        end

    end

    methods (Static, Access = protected)
        
        function vec = vectorImage_(signedPerm, vec)
        % Signed permutation of a column vector
        % equivalent to SignedPerm.matrix(signedPerm) * vec
            assert(length(signedPerm) == length(vec));
            vec(abs(signedPerm)) = vec .* sign(signedPerm(:));
        end
        
        function M = selfAdjointMatrixImage_(signedPerm, M)
        % Returns the image under the action of signedPerm on the 
        % columns and rows of M, i.e.
        % SignedPerm.matrix(perm)*M*SignedPerm.matrix(perm)'
            assert(length(signedPerm) == size(M, 1));
            assert(length(signedPerm) == size(M, 2));
            minusSign = find(signedPerm < 0);
            if length(minusSign) > 0
                M(minusSign, :) = -M(minusSign, :);
                M(:, minusSign) = -M(:, minusSign);
            end
            M(abs(signedPerm), abs(signedPerm)) = M;
        end
            
        function mat = toMatrix_(signedPerm)
        % Returns the signed permutation matrix corresponding to the given
        % signed permutation such that matrix multiplication is
        % compatible with composition of permutations, i.e. 
        % SignedPerm.matrix(SignedPerm.compose(x, y)) = 
        % SignedPerm.matrix(x) * SignedPerm.matrix(y)
            n = length(signedPerm);
            mat = sparse(abs(signedPerm), 1:n, sign(signedPerm), n, n);
        end

        function perm = toPermutation_(signedPerm)
            n = length(signedPerm);
            perm = zeros(1, 2*n);
            for i = 1:length(signedPerm)
                im = signedPerm(i);
                if im > 0
                    shift = [1 2];
                else
                    shift = [2 1];
                end
                perm((i-1)*2 + [1 2]) = (abs(im)-1)*2 + shift;
            end
        end

        function p = findMovedPoint_(signedPerm)
            for i = 1:length(signedPerm)
                if i ~= signedPerm(i)
                    p = i;
                    return
                end
            end
            p = [];
        end
        
    end
    
end
