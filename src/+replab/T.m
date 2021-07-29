function grp = T(n)
% Returns the standard torus group of given rank
%
% Args:
%   n (integer): Torus group rank
%
% Returns:
%   `+replab.TorusGroup`: Standard torus group
    grp = replab.TorusGroup(zeros(0, n));
end
