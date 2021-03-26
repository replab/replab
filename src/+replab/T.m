function grp = T(n)
% Returns the standard torus group of given rank
%
% Args:
%   n (integer): Torus group rank
%
% Returns:
%   `+replab.StandardTorusGroup`: Torus group
    grp = replab.StandardTorusGroup(n);
end
