classdef Perm
    
    methods (Static)
        
        function mat = matrix(perm)
        % Returns the permutation matrix corresponding to the given permutation
        % such that matrix multiplication is compatible with composition of
        % permutations, i.e. 
        % Perm.matrix(Perm.compose(x, y)) = Perm.matrix(x) * Perm.matrix(y)
            n = length(perm);
            mat = sparse(perm, 1:n, ones(1, n), n, n);
        end
        
        function z = compose(x, y)
        % Returns the composition of x and y
            assert(length(x) == length(y));
            z = x(y);
        end
        
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
        
        function x = identity(n)
        % Returns the identity permutation
            x = 1:n;
        end
        
        function y = inverse(x)
        % Returns the permutation y such that Perm.compose(x, y) == identity
            n = length(x);
            y = zeros(1, n);
            y(x) = 1:n;
        end
        
        function x = random(n)
        % Returns a random permutation acting on 1..n
            x = randperm(n);
        end
        
        function p = fromCycles(n, varargin)
        % Constructs a permutation from a product of cycles, each
        % cycle being a row vector, and the sequence cycles being
        % given as variable arguments
            import replab.*
            p = 1:n;
            for i = length(varargin):-1:1
                cycle = varargin{i};
                % cycle 2 3 1 means that 2 -> 3, 3 -> 1, 1 -> 2
                cycleImage = [cycle(2:end) cycle(1)];
                newEl = 1:n;
                newEl(cycle) = cycleImage;
                p = Perm.compose(newEl, p);
            end
        end
        
    end
    
end
