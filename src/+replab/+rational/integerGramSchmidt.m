function M = integerGramSchmidt(M)
% Performs the integer Gram-Schmidt orthogonalization procedure on the rows of an integer matrix
%
% Aborts the process if the magnitude of the coefficients grows bigger than an internal parameter.
%
% Args:
%   M (integer matrix): Matrix to orthogonalize, must be full row rank
%
% Returns:
%   double matrix or []: Result of the orthogonalization process
    maxabs = 60; % maximum coefficient size
    M(1,:) = M(1,:)/replab.rational.vecgcd(M(1,:));
    for i = 2:size(M, 1)
        Mi = M(i,:);
        for j = 1:i-1
            Mj = M(j,:);
            Mi = Mi*dot(Mj,Mj) - dot(Mj,Mi)*Mj;
            g = replab.rational.vecgcd(Mi);
            Mi = Mi/g;
            if max(abs(Mi)) > maxabs
                M = [];
                return
            end
        end
        M(i,:) = Mi;
    end
end
