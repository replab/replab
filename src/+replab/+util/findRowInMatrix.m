function ind = findRowInMatrix(row, matrix)
% Finds the rows of a matrix that match a given row
%
% Args:
%   row (double(1, n)): Row vector to match
%   matrix (double(\*, n)): Matrix to search
%
% Returns:
%   integer(1,\*): Row indices that match
    if isempty(matrix)
        ind = [];
        return
    end

    ind = find(matrix(:, 1) == row(1))';
    nCols = size(matrix, 2);
    for c = 2:nCols
        ind = ind(matrix(ind, c) == row(c));
        switch length(ind)
          case 0
            return
          case 1
            if any(matrix(ind, c+1:end) ~= row(c+1:end))
                ind = [];
            end
            return
          case 2
            if any(matrix(ind(2), c+1:end) ~= row(c+1:end))
                ind = ind(1);
            end
            if any(matrix(ind(1), c+1:end) ~= row(c+1:end))
                ind = ind(2:end);
            end
            return
        end
    end
end
