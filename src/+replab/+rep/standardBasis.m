function U = standardBasis(d)
% Returns an orthogonal integer basis of the Euclidean space that including the vector of all ones
%
% Args:
%   d (integer): Dimension
%
% Returns:
%   integer matrix: Matrix whose first row is the vector of all ones,
%                   with all rows orthogonal to each other
    F = factor(d);
    U = 1;
    for i = 1:length(F)
        f = F(i);
        Ui = zeros(f, f);
        Ui(1, :) = ones(f, 1);
        for j = 2:f
            Ui(j,j-1) = f-j+1;
            Ui(j,j:end) = -1;
        end
        U = kron(U, Ui);
    end
end
