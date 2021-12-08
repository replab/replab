function sub = derivedSubgroup(group)
% Computes the derived subgroup of a permutation group
%
% Args:
%   `+replab.PermutationGroup`: Permutation group
%
% Returns:
%   `+replab.PermutationGroup`: Derived subgroup
    nG = group.nGenerators;
    n = group.domainSize;
    chain = replab.bsgs.Chain(n);
    generators = {};
    for i = 1:nG
        gi = group.generator(i);
        for j = 1:nG
            gj = group.generator(j);
            cm = group.composeWithInverse(group.compose(gi, gj), group.compose(gj, gi));
            if chain.stripAndAddStrongGenerator(cm)
                generators{1, end+1} = cm;
                chain.deterministicSchreierSims;
            end
        end
    end
    % compute the normal closure
    toCheck = generators;
    while ~isempty(toCheck)
        h = toCheck{end};
        toCheck = toCheck(1:end-1);
        for i = 1:nG
            gi = group.generator(i);
            cm = group.leftConjugate(gi, h);
            if chain.stripAndAddStrongGenerator(cm)
                generators{1, end+1} = cm;
                toCheck{1, end+1} = cm;
                chain.deterministicSchreierSims;
            end
        end
    end
    chain.makeImmutable;
    sub = replab.PermutationGroup(group.domainSize, generators, 'order', chain.order, 'chain', chain);
end
