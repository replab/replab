function relators = relatorsForPermutationGroup(group, names)
    permList = @(p) ['Inverse(PermList([' strjoin(arrayfun(@num2str, p, 'uniform', 0), ', ') ']))'];
    line1 = ['G := GroupByGenerators([' strjoin(cellfun(permList, group.generators, 'uniform', 0), ', ') ']);;'];
    line2 = 'F := Image(IsomorphismFpGroupByGenerators(G, GeneratorsOfGroup(G)));;';
    line3 = 'GeneratorsOfGroup(F);';
    line4 = 'RelatorsOfFpGroup(F);';
    lines = strjoin({line1 line2 line3 line4}, '\n');
    tfile = tempname();
    fid = fopen(tfile, 'wt');
    fprintf(fid, lines);
    fclose(fid);
    [status, result] = system(['~/software/gap4/bin/gap.sh -q <' tfile]);
    delete(tfile);
    outs = strsplit(result, '\n');
    gapNames = strtrim(outs{1});
    assert(gapNames(1) == '[' && gapNames(end) == ']');
    gapNames = gapNames(2:end-1);
    gapNames = strsplit(gapNames, ',');
    gapNames = cellfun(@strtrim, gapNames, 'uniform', 0);
    relators = strtrim(outs{2});
    assert(relators(1) == '[' && relators(end) == ']');
    relators = relators(2:end-1);
    relators = strsplit(relators, ',');
    relators = cellfun(@strtrim, relators, 'uniform', 0);
    for i = 1:length(relators)
        r = relators{i};
        for j = 1:length(names)
            r = strrep(r, sprintf('F%d', j), names{j});
        end
        relators{i} = r;
    end
end
