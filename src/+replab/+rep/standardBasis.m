function U = standardBasis(d)
% Returns an orthogonal integer basis of R^d including the vector of all ones
%
% Basis elements are columns of U
    F = factor(d);
    U = 1;
    for i = 1:length(F)
        f = F(i);
        Ui = zeros(f, f);
        Ui(:, 1) = ones(f, 1);
        for j = 2:f
            Ui(j-1,j) = f-j+1;
            Ui(j:end,j) = -1;
        end
        U = kron(U, Ui);
    end
    % new convention with basis vectors in rows
    U = U';
end
