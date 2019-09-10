function [l r] = brackets(obj)
% Returns the pair {} for cell arrays or [] for other objects
    if iscell(obj)
        l = '{'; r = '}';
    else
        l = '['; r = ']';
    end
end
