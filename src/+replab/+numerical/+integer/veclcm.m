function l = veclcm(v)
% Computes the lcm of all elements in a vector
%
% Args:
%   v (integer(1,\*)): Integer vector
%
% Returns:
%   integer: Least common multiple
    if isempty(v)
         l = 1;
    else
        v = unique(v);
        l = v(1);
        for i = 2:length(v)
            l = lcm(l, v(i));
        end
    end
end