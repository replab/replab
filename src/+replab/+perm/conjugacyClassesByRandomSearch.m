function classes = conjugacyClassesByRandomSearch(group)
% Generates the conjugacy classes of a group
%
% Based on the randomized method

% Args:
%   group (`+replab.PermutationGroup`): Permutation group to compute the conjugacy classes of
%
% Returns:
%   classImages (cell(1,\*) of `+replab.ConjugacyClass`): All conjugacy classes
    classes = {};
    remains = group.order;
    if group.order > 2^52
        factor = group.order/2^52;
    else
        factor = 1;
    end
    bar = replab.infra.repl.ProgressBar(group.order);
    while remains > 0
        g = group.sample;
        gn = g;
        n = 1;
        eo = double(group.elementOrder(g));
        while n < eo
            bar.step(group.order - remains);
            newClass = true;
            for i = 1:length(classes)
                c = classes{i};
                r = c.representative;
                if isequal(replab.Permutation.cycleStructure(r), replab.Permutation.cycleStructure(gn))
                    if c.contains(g)
                        newClass = false;
                        break
                    end
                end
            end
            if newClass
                c = replab.ConjugacyClass(group, gn);
                classes{1,end+1} = c;
                remains = remains - c.size;
            end
            gn = gn(g);
            n = n + 1;
        end
    end
    bar.finish;
end
