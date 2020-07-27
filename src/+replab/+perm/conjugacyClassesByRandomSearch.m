function classes = conjugacyClassesByRandomSearch(group)
% Generates the conjugacy classes of a group
%
% Based on the randomized method, see p. 215 of
% A. Seress, Permutation Group Algorithms. Cambridge University Press, 2003.
%
% and
%
% M. Jerrum, "Computational Pólya theory," in Surveys in combinatorics, 1995, 1995, pp. 103–118.
% Args:
%   group (`+replab.PermutationGroup`): Permutation group to compute the conjugacy classes of
%
% Returns:
%   classImages (cell(1,\*) of `+replab.ConjugacyClass`): All conjugacy classes
    if group.order > 2^52
        factor = group.order/2^52;
    else
        factor = 1;
    end
    bar = replab.infra.repl.ProgressBar(group.order);
    classes = {};
    remains = group.order;
    g = group.sample;
    while remains > 0
        bar.step(group.order - remains);
        newClass = true;
        for i = 1:length(classes)
            c = classes{i};
            r = c.representative;
            if isequal(replab.Permutation.cycleStructure(r), replab.Permutation.cycleStructure(g))
                if c.contains(g)
                    newClass = false;
                    break
                end
            end
        end
        if newClass
            c = replab.ConjugacyClass.make(group, g);
            classes{1,end+1} = c;
            remains = remains - c.size;
        end
        g = c.representativeCentralizer.sample;
    end
    assert(remains == 0);
    bar.finish;
end
