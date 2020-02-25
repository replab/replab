function res = areAllEqual(vec)
% Returns whether all the elements in a vector are equal under isequal
    res = true;
    for i = 2:length(vec)
        if ~isequal(vec{1}, vec{i})
            res = false;
            return
        end
    end
end
