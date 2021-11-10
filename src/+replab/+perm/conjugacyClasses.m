function res = conjugacyClasses(group)
% Computes the conjugacy classes of a permutation group
%
% Args:
%   group (`+replab.PermutationGroup`): Permutation group
%
% Returns:
%   `+replab.ConjugacyClasses`: Conjugacy classes
    if group.order < 100000
        C = replab.perm.conjugacyClassesByOrbits(group);
        n = length(C);
        classes = cell(1, n);
        for i = 1:n
            cl = sortrows(C{i}');
            classes{i} = replab.ConjugacyClass(group, cl(1,:));
        end
    else
        classes = replab.perm.conjugacyClassesByRandomSearch(group);
    end
    reps = zeros(length(classes), group.domainSize);
    for i = 1:length(classes)
        reps(i,:) = classes{i}.representative;
    end
    [~, I] = sortrows(reps);
    res = replab.ConjugacyClasses(self, classes(I)); % sort by minimal representative
end
