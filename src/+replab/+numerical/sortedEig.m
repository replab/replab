function [U D] = sortedEig(A, direction, useAbs)
% Returns the eigenvalues, sorted by magnitude
    if nargin < 2
        direction = 'ascend';
        useAbs = true;
    end
    [U D] = eig(A); % eigenvalue decomposition
    D = diag(D);
    if useAbs
        [~, I] = sort(abs(D), direction);
    else
        [~, I] = sort(D, direction);
    end
    D = diag(D(I)); % reorder accordingly
    U = U(:, I);
end
