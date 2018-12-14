classdef Perm
% Helper methods for permutations
    
    methods (Static) % CONSTRUCTION METHODS
        
        function x = identity(n)
        % Returns the identity permutation
            x = 1:n;
        end
        
        function p = fromCycles(n, varargin)
        % Constructs a permutation from a product of cycles, each
        % cycle being a row vector, and the sequence cycles being
        % given as variable arguments
            p = 1:n;
            for i = length(varargin):-1:1
                cycle = varargin{i};
                % cycle 2 3 1 means that 2 -> 3, 3 -> 1, 1 -> 2
                cycleImage = [cycle(2:end) cycle(1)];
                newEl = 1:n;
                newEl(cycle) = cycleImage;
                p = replab.Perm.compose(newEl, p);
            end
        end
        
        function x = random(n)
        % Returns a random permutation acting on 1..n
            x = randperm(n);
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
            z = x(y);
        end
        
        function y = pow(x, e)
        % Computes y = x^e by repeated squaring
            n = length(x);
            if e < 0
                y = replab.Perm.pow(replab.Perm.inverse(x), -e);
            elseif e == 0
                y = replab.Perm.identity(n);
            else
                y = replab.Perm.identity(n);
                while e > 1
                    if mod(e, 2) == 0 % n even
                        x = replab.Perm.compose(x, x);
                        e = e / 2;
                    else
                        y = replab.Perm.compose(x, y);
                        x = replab.Perm.compose(x, x);
                        e = (e - 1)/2;
                    end
                end
                y = replab.Perm.compose(x, y);
            end
        end
        
        function y = inverse(x)
        % Returns the permutation y such that Perm.compose(x, y) == identity
            n = length(x);
            y = zeros(1, n);
            y(x) = 1:n;
        end
        
    end

    methods (Static) % REPRESENTATIONS
        
        function mat = matrix(perm)
        % Returns the permutation matrix corresponding to the given permutation
        % such that matrix multiplication is compatible with composition of
        % permutations, i.e. 
        % Perm.matrix(Perm.compose(x, y)) = Perm.matrix(x) * Perm.matrix(y)
            n = length(perm);
            mat = sparse(perm, 1:n, ones(1, n), n, n);
        end
        
    end
        
    methods (Static) % ACTIONS
        
        function i = image(perm, p)
        % Returns the image i of the point p under the permutation perm
            i = perm(p);
        end
        
        function vec = vectorAction(perm, vec)
        % Permutation of a column vector
        % equivalent to Perm.matrix(perm) * vec
            assert(length(perm) == length(vec));
            vec(perm) = vec;
        end
        
        function M = selfAdjointMatrixAction(perm, M)
        % Returns the image under the action of perm on the columns and rows of M
        % i.e. Perm.matrix(perm)*M*Perm.matrix(perm)'
            n = length(perm);
            assert(size(M, 1) == n);
            assert(size(M, 2) == n);
            M(perm, perm) = M;
        end
        
    end
    
end
