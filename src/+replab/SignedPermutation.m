classdef SignedPermutation
% This class regroups functions that are useful when dealing with signed permutations.
    
    methods (Static)

        function perm = toPermutation(signedPerm)
        % Returns the permutation corresponding to the given signed permutation where the permutation acts on the list [-d,..,-1,1,..,d]
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
        %
        % See `toPermutation`
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
        % Returns the signed permutation matrix corresponding to the given signed permutation
        %
        % Such that matrix multiplication is
        % compatible with composition of permutations, i.e.
        % S.toMatrix(S.compose(x, y)) =
        % S.toMatrix(x) * S.toMatrix(y)
        % where S = SignedPermutations(domainSize)
            mat = full(replab.SignedPermutation.toSparseMatrix(signedPerm));
        end

        function b = isSignedPermutationMatrix(mat)
        % Returns true when "mat" is a signed permutation matrix, i.e. a monomial matrix with nonzero entries equal to +1 or -1
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

        function mat = toSparseMatrix(signedPerm)
            n = length(signedPerm);
            mat = sparse(abs(signedPerm), 1:n, sign(signedPerm), n, n);
        end


        function signedPerm = fromMatrix(mat)
        % Returns the signed permutation corresponding to the given matrix representation or throws an error
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

end
