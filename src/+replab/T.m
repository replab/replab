function grp = T(n)
% Returns the torus group acting of a given dimension
%
% Args:
%   n (integer): Torus group dimension
%
% Returns:
%   `+replab.TorusGroup`: Torus group
    grp = replab.TorusGroup(n);
end
