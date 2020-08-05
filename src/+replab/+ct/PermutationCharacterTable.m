function ct = PermutationCharacterTable(group)
% Generates the character table for a permutation group with integer characters
% 
% Args:
%   group (`+replab.PermutationGroup`): permutation group
%
% Returns:
%   ct (`+replab.CharacterTable`)
    assert(isa(group, 'replab.PermutationGroup'));
    ord = double(group.order);
    decomp = group.naturalRep.decomposition.nice;
    irreps = decomp.components;
    classes = group.conjugacyClasses;
    k = length(classes);
    for i = 1:k
        % verify that g^n is conjugate to g for all g \in G and n \in Z with n and (g) coprime
        % https://math.stackexchange.com/questions/2792741/classification-of-groups-with-integer-valued-characters
        cl = classes{i};
        g = cl.representative;
        eo = group.elementOrder(g);
        for j = 2:eo-1
            if gcd(j, eo) == 1
                if ~cl.contains(group.composeN(g, j))
                    error('Group does not have integer characters');
                end
            end
        end
    end
    ccreps = cell(1, length(classes));
    cclens = cell(1, length(classes));
    for i = 1:k
        ccreps{i} = classes{i}.representative;
        cclens{i} = double(classes{i}.nElements);
    end
    nirreps = length(irreps);
    chars = cell(nirreps, k);
    for i = 1:nirreps
        irrep = irreps{i};
        chars(i,:) = cellfun(@(x) trace(irrep.image(x)), ccreps, 'UniformOutput', false);
    end
    chars = cell2mat(chars);
    start_irreps = 1;
    while nirreps ~= k
        for i = start_irreps:nirreps
            for j = start_irreps:nirreps
                new_rep = irreps{i}.kron(irreps{j});
                ss = 0;
                for n = 1:k
                    char = trace(new_rep.image(ccreps{n}));
                    ss = ss + cclens{n}*(char)^2;
                end
                if ss ~= ord
                    new_irreps = new_rep.decomposition.nice;
                    start_irreps = nirreps + 1;
                    for m = 1:new_irreps.nComponents
                        new_irrep = new_irreps.component(m);
                        new_char = cellfun(@(x) trace(new_irrep.image(x)), ...
                                            ccreps, 'UniformOutput', false);
                        new_char = cell2mat(new_char);
                        if ~any(all(abs(chars - new_char) < 1/2, 2))
                            irreps{nirreps + 1} = new_irrep;
                            chars(nirreps + 1, :) = new_char;
                            nirreps = nirreps + 1;
                        end
                    end
                end
            end
        end
    end
    chars = num2cell(round(chars));
    chars = cellfun(@num2str, chars, 'UniformOutput', false);
    gens = group.generators;
    irrepExp = cell(1, nirreps);
    for i = 1:nirreps
        genImages = cellfun(@(x) irreps{i}.image(x), gens, 'UniformOutput', false);
        genImages = cellfun(@(x) cellfun(@num2str, num2cell(x), 'UniformOutput', false), genImages, 'UniformOutput', false);
        irrepExp{i} = genImages;
    end
    ct = replab.CharacterTable.make(group, classes, [], chars, [], irrepExp);
end

