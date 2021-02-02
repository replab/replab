function [V,D,W] = realeig(A)
% Eigendecomposition of real matrices using real matrices
%
% The matrices returned obey `` A*V = V*D `` and `` W'*A = D*W' ``.
%
% This is similar to the use of ``eig`` followed by ``cdf2rdf``, except that we perform additional
% post-processing.
%
% * We make sure that 2x2 complex blocks in ``D`` have the form ``[a -b; b a]`` with ``b > 0``.
% * We sort the eigenvalues in increasing order of their real parts.
%
% Args:
%   A (double(d,d)): Square real matrix, not necessarily symmetric
%
% Returns
% -------
%   V: double(d,d)
%     Right eigenvectors
%   D: double(d,d), sparse
%     Eigenvalues, with complex eigenvalues encoded in 2x2 real blocks
%   W: double(d,d)
%     Left eigenvectors (in column form)
    assert(isreal(A));
    d = size(A, 1);
    [V,D,W] = eig(A,'vector');
    D = reshape(D, 1, []);
    complexInd = reshape(find(imag(D) ~= 0), 1, []);
    realInd = reshape(find(imag(D) == 0), 1, []);
    blocks = [realInd; zeros(1, length(realInd))];
    realPart = D(realInd);
    D = sparse(1:d, 1:d, D, d, d); % = diag(D)
    for j = complexInd(1:2:end)
        % our block has D of the form
        % [a+i*b      ]             [d1   ]
        % [      a-i*b]  which is   [   d2]
        % and will be transformed into (as in the convention of cdf2rdf)
        % [ a b
        %  -b a]
        % thus we want b < 0 because that is our convention (as in Pauli's sigma_y has a negative sign up-right)
        d1 = D(j,j);
        d2 = D(j+1,j+1);
        assert(d1 == conj(d2));
        a = real(d1);
        b = imag(d1);
        if b > 0
            V(:,[j j+1]) = V(:,[j+1 j]);
            W(:,[j j+1]) = W(:,[j+1 j]);
            b = -b;
        end
        blocks = [blocks [j;j+1]];
        realPart = [realPart a];
        V(:,j:j+1) = V(:,j:j+1) * [1 -1i; 1 1i]/2;
        W(:,j:j+1) = W(:,j:j+1) * [1 -1i; 1 1i]/2;
        D(j:j+1,j:j+1) = [ a b
                          -b a];
    end
    assert(isreal(V) && isreal(W));
    [~, I] = sort(realPart);
    blocks = reshape(blocks(:, I), 1, []);
    blocks = blocks(blocks > 0);
    V = V(:, blocks);
    W = W(:, blocks);
    D = D(blocks, blocks);
end
