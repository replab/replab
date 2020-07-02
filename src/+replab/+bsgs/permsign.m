function sign = permsign(perm)
% Returns the sign of a given permutation
%
% Args:
%   perm (permutation): Vector representing a permutation (e.g. [3 2 1 4])
%
% Returns:
%   integer: Sign of the permutation
    x = perm;
    n = length(x);
    oddOrEven = 0; %Records whether the total permutation is odd or even
    for i = 1:n
        if x(i) == 0 || x(i) == i %Skip over one cycles and numbers that have been cycled through
            continue
        end
        cycleSize = -1; %The first element in a cycle isn't counted
        j = i;
        while x(j) ~= 0
            pHold = x(j);
            x(j) = 0;
            j = pHold;
            cycleSize = cycleSize + 1;
        end
        if cycleSize > 0
            oddOrEven = oddOrEven + cycleSize; %At the end, this will match the parity (even/odd) of the permuation
        end
    end
    sign = (-1)^mod(round(oddOrEven),2); %Sign of permutation
end
