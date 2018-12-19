classdef SignedPermutations < replab.cat.Group
   
    properties (SetAccess = protected)
        n;
        canEqv = true;
        canHash = true;
        canSample = true;
    end
    
    methods
        
        function self = SignedPermutations(n)
            self.n = n;
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
                                         @(g) replab.SignedPermutations.findMovedPoint(g));
        end
        
        function A = vectorAction(self, field)
            n = self.n;
            desc = sprintf('Signed permutations acting on %d dimensional vectors in %s', n, field);
            P = replab.Vectors(n, field);
            A = replab.cat.ActionFun(desc, self, P, ...
                                     @(g, p) replab.SignedPermutations.vectorImage(g, p));
        end
        
        function A = selfAdjointMatrixAction(self, field)
            n = self.n;
            desc = sprintf('Signed permutations acting on %d x %d self-adjoint matrices in %s', n, field);
            P = replab.SelfAdjointMatrices(n, field);
            A = replab.cat.ActionFun(desc, self, P, ...
                                     @(g, p) replab.SignedPermutations.selfAdjointMatrixImage(g, p));
        end
        
        function phi = permutationIsomorphism(self)
            fun = @(x) replab.SignedPermutations.toPermutation(x);
            target = replab.Permutations(2*self.n);
            phi = replab.cat.GroupMorphismFun(fun, self, target);
        end

    end
    
    methods
       
        function law_check_toMatrix_fromMatrix_D(self, signedPerm)
            M = replab.SignedPermutations.toMatrix(signedPerm);
            P = replab.SignedPermutations.fromMatrix(M);
            self.assertEqv(signedPerm, P);
        end
        
    end

    methods (Static)
        
        function vec = vectorImage(signedPerm, vec)
        % Signed permutation of a column vector
        % equivalent to SignedPerm.matrix(signedPerm) * vec
            assert(length(signedPerm) == length(vec));
            vec(abs(signedPerm)) = vec .* sign(signedPerm(:));
        end
        
        function M = selfAdjointMatrixImage(signedPerm, M)
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
            
        function mat = toMatrix(signedPerm)
        % Returns the signed permutation matrix corresponding to the given
        % signed permutation such that matrix multiplication is
        % compatible with composition of permutations, i.e. 
        % SignedPerm.matrix(SignedPerm.compose(x, y)) = 
        % SignedPerm.matrix(x) * SignedPerm.matrix(y)
            n = length(signedPerm);
            mat = sparse(abs(signedPerm), 1:n, sign(signedPerm), n, n);
        end
        
        function signedPerm = fromMatrix(mat)
            signedPerm = [];
            n = size(mat, 1);
            [I J S] = find(mat);
            if length(I) ~= n
                return
            end
            I = I';
            J = J';
            S = S';
            sI = sort(I);
            [sJ IJ] = sort(J);
            if ~isequal(sI, 1:n) || ~isequal(sJ, 1:n)
                return
            end
            signedPerm = I.*S;
            signedPerm = signedPerm(IJ);
        end

        function perm = toPermutation(signedPerm)
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

        function p = findMovedPoint(signedPerm)
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
