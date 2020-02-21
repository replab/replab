function irreps = splitAll(rep, context)
% Decomposes fully the given representation into subrepresentations
    irreps = cell(1, 0);
    todo = {rep};
    while ~isempty(todo)
        h = todo{1};
        if isequal(h.isIrreducible, true)
            % head of list is irreducible, remove it
            irreps{1,end+1} = h;
            todo = todo(2:end);
        else
            res = replab.irreducible.split(h, context);
            if length(res) == 1
                todo{1} = res;
            else
                todo = horzcat(todo(2:end), res);
            end
        end
    end
end
