classdef SignedPerm
% Helper methods for signed permutations

    methods (Static) % CONSTRUCTION METHODS
               
        function x = identity(n)
        % Returns the identity signed permutation
            x = 1:n;
        end
                
        function x = random(n)
        % Returns a random permutation acting on 1..n
            s = randi([0 1], 1, n)*2-1;
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
        
        function perm = toPerm(signedPerm)
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

        function mat = toMatrix(signedPerm)
        % Returns the signed permutation matrix corresponding to the given
        % signed permutation such that matrix multiplication is
        % compatible with composition of permutations, i.e. 
        % SignedPerm.matrix(SignedPerm.compose(x, y)) = 
        % SignedPerm.matrix(x) * SignedPerm.matrix(y)
            n = length(signedPerm);
            mat = sparse(abs(signedPerm), 1:n, sign(signedPerm), n, n);
        end
        
    end
        
    methods (Static) % ACTIONS
        
        function i = image(signedPerm, p)
        % Returns the image i of the point p under the signed permutation perm
            i = signedPerm(abs(p))*sign(p);
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

        function i1 = invImage(i, signedPerm)
            i1 = find(abs(signedPerm) == abs(i));
            i1 = i1 * sign(signedPerm(i1)) * sign(i);
        end
        
        function vec = vectorImage(signedPerm, vec)
        % Signed permutation of a column vector
        % equivalent to SignedPerm.matrix(signedPerm) * vec
            assert(length(signedPerm) == length(vec));
            vec(abs(signedPerm)) = vec .* sign(perm(:));
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
            M(abs(gp), abs(gp)) = M;
        end
        
    end
    
end
