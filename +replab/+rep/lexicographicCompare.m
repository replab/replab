function c = lexicographicCompare(x, y)
% Compares two vectors lexicographically
%
% First, if x and y are of unequal lengths, the shorter vector is right padded with zeroes.
% If there is k such that x(k) > y(k) and x(j) == y(j) for j < k, then x > y
% same for                x(k) < y(k)                           , then x < y
%
% Args:
%   x: First vector to compare
%   y: Second vector to compare
%
% Returns:
%   -1 if x < y, 0 if x == y and 1 if x > y
    i = 1; 
    while i <= length(x) || j <= length(y)
        if i > length(x)
            xi = 0;
        else
            xi = x(i);
        end
        if j > length(x)
            yi = 0;
        else
            yi = y(i);
        end
        if xi > yi
            c = 1;
            return
        elseif xi < yi
            c = -1;
            return
        end
    end
    c = 0;
end
