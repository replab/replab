classdef Permutation
% This class regroups functions that are useful when dealing with permutations.
%
% In RepLAB, permutations are represented using row vectors of integers, those integers being
% encoded using the type ``double``; so this class does not have any instances.

    methods (Static)

        function s = sign(perm)
        % Returns the sign of a given permutation
        %
        % Args:
        %   perm (permutation): Vector representing a permutation (e.g. [3 2 1 4])
        %
        % Returns:
        %   integer: Sign of the permutation
            x = perm;
            n = length(x);
            oddOrEven = 0; %Records whether the total permutation is odd or even
            for i = 1:n
                if x(i) == 0 || x(i) == i %Skip over one cycles and numbers that have been cycled through
                    continue
                end
                cycleSize = -1; %The first element in a cycle isn't counted
                j = i;
                while x(j) ~= 0
                    pHold = x(j);
                    x(j) = 0;
                    j = pHold;
                    cycleSize = cycleSize + 1;
                end
                if cycleSize > 0
                    oddOrEven = oddOrEven + cycleSize; %At the end, this will match the parity (even/odd) of the permuation
                end
            end
            s = (-1)^mod(round(oddOrEven),2); %Sign of permutation
        end

        function p = sorting(array, greaterFun)
        % Returns the permutation that sorts a cell array using a custom comparison function
        %
        % Args:
        %   array: A (row or column) vector array containing data
        %          of arbitrary type
        %   greaterFun: A function handle that compares elements such that
        %               ``greaterFun(x, y) == true`` when x > y
        %               and false otherwise
        % Returns:
        %   A permutation ``p`` such that ``sorted = array(p)``
            n = length(array);
            p = 1:n;
            inc = round(n/2);
            while inc > 0
                for i = (inc+1):n
                    tmp = p(i);
                    j = i;
                    if isa(array, 'cell')
                        while (j >= inc+1) && greaterFun(array{p(j-inc)}, array{tmp})
                            p(j) = p(j-inc);
                            j = j - inc;
                        end
                    else
                        while (j >= inc+1) && greaterFun(array(p(j-inc)), array(tmp))
                            p(j) = p(j-inc);
                            j = j - inc;
                        end
                    end
                    p(j) = tmp;
                end
                if inc == 2
                    inc = 1;
                else
                    inc = round(inc/2.2);
                end
            end
        end

        function mat = toSparseMatrix(perm)
        % Returns the sparse permutation matrix corresponding to the given permutation
        %
        % The returned matrix is such that matrix multiplication is compatible with composition of
        % permutations, i.e. for ``P = replab.SymmetricGroup(domainSize)`` we have
        % ``P.toMatrix(P.compose(x, y)) = P.toMatrix(x) * P.toMatrix(y)``
        %
        % Args:
        %   perm (permutation): Permutation
        %
        % Returns:
        %   The sparse permutation matrix corresponding to ``perm``.
            n = length(perm);
            mat = sparse(perm, 1:n, ones(1, n), n, n);
        end

        function mat = toMatrix(perm)
        % Returns the permutation matrix corresponding to the given permutation
        %
        % The returned matrix is such that matrix multiplication is compatible with composition of
        % permutations, i.e. for ``P = replab.Permutations(domainSize)`` we have
        % ``P.toMatrix(P.compose(x, y)) = P.toMatrix(x) * P.toMatrix(y)``
        %
        % Args:
        %   perm (permutation): Permutation
        %
        % Returns:
        %   The permutation matrix corresponding to ``perm``.
            mat = full(replab.Permutation.toSparseMatrix(perm));
        end

        function perm = fromMatrix(mat)
        % Returns the permutation corresponding to the given matrix representation
        %
        % See `+replab.Permutation.toMatrix`
        %
        % Args:
        %   mat: A permutation matrix.
        %
        % Returns:
        %   The permutation corresponding to matrix ``mat``.
        %
        % Raises:
        %   Error: if ``mat`` is not a permutation matrix, throws an error
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
