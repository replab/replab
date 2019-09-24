function [intBasisOpt isOrtho] = integerRowSpan(basis)
% Attempts integer basis recovery from the given unitary basis
%
% Tries to find a rational matrix that has the same row span as the given basis `basis`. 
% It uses the `rat` Matlab function that uses truncated continued fraction expansions. 
% Additionally, the function attemps to perform Gram Schmidt orthogonalization over the integers
% (but will skip that step if it would make the coefficients to grow too big).
%
% The method fails when the recovered approximation has integer coefficients that are too big.
%
% It also uses the heuristic that the orthogonal projector on the original space has rational coefficients.
%
% Args:
%   basis (double matrix): Basis of dimension nRows x nCols with nRows <= nCols
%
% Returns
% -------
%   intBasisOpt: integer matrix or []
%     Integer decomposition of the basis, or [] if the recovery failed
%   isOrtho: logical
%     Whether the integer Gram-Schmidt process succeeded and the basis vectors
%     are actually orthogonal
%
% Raises:
%   An error if the basis is complex
    assert(isreal(basis), 'The basis must be real');
    nRows = size(basis, 1);
    Q = basis*basis';
    % make the given basis orthonormal if that is not the case
    if replab.isNonZeroMatrix(Q - eye(nRows), replab.Settings.doubleEigTol)
        if isdiag(Q)
            % special case for diagonal correction, use square root
            T = diag(1./sqrt(diag(Q)));
        else
            T = inv(chol(Q, 'lower'));
        end
        basis = T * basis;
    end
    % compute the projector
    P = basis'*basis;
    % find a well conditioned basis, see
    % https://www.mathworks.com/matlabcentral/answers/108835-how-to-get-only-linearly-independent-rows-in-a-matrix-or-to-remove-linear-dependency-b-w-rows-in-a-m
    [~,~,jb] = qr(P, 'vector');
    jb = jb(1:nRows);
    U = P(jb, :);
    assert(length(jb) == nRows);
    [num den] = replab.rational.attemptRecoverRational(U);
    if isempty(num)
        intBasisOpt = [];
        isOrtho = [];
    else
        afterGS = replab.rational.integerGramSchmidt(num);
        if isempty(afterGS)
            intBasisOpt = num;
            isOrtho = false;
        else
            intBasisOpt = afterGS;
            isOrtho = true;
        end
    end
end
