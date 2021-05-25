function K = computeKronecker(ct)
% Computes the generalized Kronecker coefficients for the given character table
%
% Args:
%   ct (`+replab.CharacterTable`): Real or complex character table

% Returns:
%   integer(\*,\*,\*): Kronecker coefficients
    n = ct.classes.nClasses;
    K = zeros(n,n,n);
    for j = 1:n
        cj = ct.character(j);
        for k = j:n
            ck = ct.character(k);
            K(:,j,k) = ct.multiplicities(cj*ck);
            K(:,k,j) = K(:,j,k);
        end
    end
end
