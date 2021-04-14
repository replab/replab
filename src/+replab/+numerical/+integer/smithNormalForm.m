function [U,S,V] = smithNormalForm(A)
% Computes the Smith normal form of an integer matrix
%
% ``[U,S,V] = replab.numerical.integer.smithNormalForm(A)`` returns integer matrices U, S, and V such that
% ``A = U*S*V'``. ``S`` is diagonal and nonnegative, ``S(i,i)`` divides ``S(i+1,i+1)`` for all i,
% and ``abs(det(U)) = 1``, ``abs(det(V)) = 1``.
%
% ``s = smith(A)`` just returns ``diag(S)``.
%
% This function is in some ways analogous to SVD.
%
% Does not pivot! so -> NaN
%
% Original code by:
% John Gilbert, 415-812-4487, December 1993
% gilbert@parc.xerox.com
% Xerox Palo Alto Research Center
%
% See https://www.researchgate.net/post/Is-there-any-computer-algorithm-that-finds-Smith-Normal-Form-of-a-given-matrix
%
% Args:
%   A (integer(m,n)): Integer matrix
%
% Returns
% -------
%   U: integer(m,m)
%     Unimodular matrix
%   S: integer(m,n)
%     Diagonal matrix
%   V: integer(n,n)
%     Unimodular matrix

    if round(A) ~= A
        error('Requires integer input.');
    end

    % This looks much like an SVD algorithm that first bidiagonalizes
    % A by Givens rotations and then chases zeros, except for
    % the construction of the 2 by 2 elementary transformation.

    [m,n] = size(A);
    S = A;
    U = eye(m);
    V = eye(n);

    % Bidiagonalize S with elementary Hermite transforms.

    for j = 1:min(m,n)
        % Zero column j below the diagonal.
        for i = j+1:m
            if S(i,j)
                % Construct an elementary Hermite transformation E
                % to zero S(i,j) by combining rows i and j.
                E = replab.numerical.integer.elementaryHermite(S(j,j),S(i,j));
                % Apply the transform to S and U.
                S([j i],:) = E * S([j i],:);
                U(:,[j i]) = U(:,[j i]) / E;
            end
        end
        % Zero row j after the superdiagonal.
        for i = j+2:n
            if S(j,i)
                % Construct an elementary Hermite transformation E
                % to zero S(j,i) by combining columns j+1 and i.
                E = replab.numerical.integer.elementaryHermite(S(j,j+1),S(j,i));
                % Apply the transform to S and V.
                S(:,[j+1 i]) = S(:,[j+1 i]) * E';
                V(:,[j+1 i]) = V(:,[j+1 i]) / E;
            end
        end
    end

    % Now S is upper bidiagonal.
    % Chase the superdiagonal nonzeros away.

    if size(S, 2) == 1
        D = [];
    elseif size(S, 1) == 1
        D = S(1, 2);
    else
        D = diag(S, 1);
    end

    while any(D)
        b = find(D, 1);
        S
        b
        % Start chasing bulge at first nonzero superdiagonal element.

        % To guarantee reduction in S(b,b), first make S(b,b) positive
        % and make S(b,b+1) nonnegative and less than S(b,b).
        if S(b,b) < 0
            S(b,:) = -S(b,:);
            U(:,b) = -U(:,b);
        end
        q = floor(S(b,b+1)/S(b,b));
        E = [1 0 ; -q 1];
        S(:,[b b+1]) = S(:,[b b+1]) * E';
        V(:,[b b+1]) = V(:,[b b+1]) / E;

        if S(b,b+1)
            % Zero the first nonzero superdiagonal element
            % using columns b and b+1, to start the bulge at S(b+1,b).
            E = replab.numerical.integer.elementaryHermite(S(b,b),S(b,b+1));
            S(:,[b b+1]) = S(:,[b b+1]) * E';
            V(:,[b b+1]) = V(:,[b b+1]) / E;
            for j = 1:min(m,n)
                if j+1 <= m
                    % Zero S(j+1,j) using rows j and j+1.
                    E = replab.numerical.integer.elementaryHermite(S(j,j),S(j+1,j));
                    S([j j+1],:) = E * S([j j+1],:);
                    U(:,[j j+1]) = U(:,[j j+1]) / E;
                end
                if j+2 <= n
                    % Zero S(j,j+2) using columns j+1 and j+2.
                    E = replab.numerical.integer.elementaryHermite(S(j,j+1),S(j,j+2));
                    S(:,[j+1 j+2]) = S(:,[j+1 j+2]) * E';
                    V(:,[j+1 j+2]) = V(:,[j+1 j+2]) / E;
                end
            end
        end
        if size(S, 2) == 1 % diag(S, 1) in a robust way
            D = [];
        elseif size(S, 1) == 1
            D = S(1, 2);
        else
            D = diag(S, 1);
        end
    end

    % Now S is diagonal.  Make it nonnegative.

    for j = 1:min(m,n)
        if S(j,j) < 0
            S(j,:) = -S(j,:);
            U(:,j) = -U(:,j);
        end
    end

    % Squeeze factors to lower right to enforce divisibility condition.

    for i = 1 : min(m,n)
        for j = i+1 : min(m,n)
            % Replace S(i,i), S(j,j) by their gcd and lcm respectively.
            a = S(i,i);
            b = S(j,j);
            [g,c,d] = gcd(a,b);
            E = [ 1 d ; -b/g a*c/g];
            F = [ c 1 ; -b*d/g a/g];
            S([i j],[i j]) = E * S([i j],[i j]) * F';
            U(:,[i j]) = U(:,[i j]) / E;
            V(:,[i j]) = V(:,[i j]) / F;
        end
    end

    U = round(U);
    V = round(V);

    if nargout <= 1
        U = diag(S);
    end

end
