function names = defaultGeneratorNames(n)
% Returns default generator names
%
% Args:
%   n (integer): Number of generators
%
% Returns:
%   cell(1,\*) of charstring: A list of generator names starting with ``x1``, ``x2``, ...
    names = arrayfun(@(i) ['x' num2str(i)], 1:n, 'uniform', 0);
end
