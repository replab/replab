classdef SignedPerm
% Helper methods for signed permutations
    
    methods (Static) % CONSTRUCTION METHODS
        
        function x = identity(n)
        % Returns the identity signed permutation
            x = 1:n;
        end
                
        function x = random(n)
        % Returns a random permutation acting on 1..n
            s = (randi(2, 1, n)-1)*2-1;
            x = randperm(n).*s;
        end
        
    end
    
    methods (Static) % PROPERTIES
        
        function b = isIdentity(perm)
        % Tests whether the given permutation is the identity
            n = length(perm);
            b = true;
            for i = 1:n
                if perm(i) ~= i
                    b = false;
                    return
                end
            end
        end

    end
    
    methods (Static) % TRANSFORMATIONS
        
        function z = compose(x, y)
        % Returns the composition of x and y
            assert(length(x) == length(y));
            z = x(abs(y)).*sign(y);
        end
        
        function y = inverse(x)
        % Returns the permutation y such that Perm.compose(x, y) == identity
            n = length(x);
            y = zeros(1, n);
            xAbs = abs(x);
            y(xAbs) = 1:n;
            invFlip = xAbs(x < 0);
            y(invFlip) = -y(invFlip);
        end
        
    end

    methods (Static) % REPRESENTATIONS
        
        function mat = matrix(signedPerm)
        % Returns the permutation matrix corresponding to the given
        % signed permutation such that matrix multiplication is
        % compatible with composition of permutations, i.e. 
        % SignedPerm.matrix(SignedPerm.compose(x, y)) = 
        % SignedPerm.matrix(x) * SignedPerm.matrix(y)
            n = length(signedPerm);
            mat = sparse(signedPerm, 1:n, sign(signedPerm), n, n);
        end
        
    end
        
    methods (Static) % ACTIONS
        
        function i = image(signedPerm, p)
        % Returns the image i of the point p under the signed permutation perm
            i = signedPerm(abs(p))*sign(p);
        end
        
        function vec = vectorAction(signedPerm, vec)
        % Signed permutation of a column vector
        % equivalent to SignedPerm.matrix(signedPerm) * vec
            assert(length(signedPerm) == length(vec));
            vec(abs(signedPerm)) = vec .* sign(perm(:));
        end
        
        function M = selfAdjointMatrixAction(signedPerm, M)
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
            M(abs(gp), abs(gp)) = M;
        end
        
    end
    
end
