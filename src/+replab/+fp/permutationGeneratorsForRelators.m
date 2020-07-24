function gens = permutationGeneratorsForRelators(names, relators)
% Computes a permutation realization of a finite group given by a presentation
%
% Note: calls the GAP system internally
%
% Args:
%   names (cell(1,\*) of charstring): Generator names
%   relators (cell(1,\*) of charstring): Relators given as explicit words
%
% Returns:
%   cell(1,\*) of permutation: Generators of a permutation group realizing the presentation
    gapRelators = cell(1, length(relators));
    gapNames = arrayfun(@(i) sprintf('f.%d', i), 1:length(names), 'uniform', 0);
    for i = 1:length(relators)
        [ok, tokens] = replab.fp.Parser.lex(relators{i}, names);
        assert(ok);
        [pos, letters] = replab.fp.Parser.word(tokens, 1);
        assert(tokens(1, pos) == replab.fp.Parser.types.END);
        gapRelators{i} = replab.fp.printLetters(letters, gapNames, '*');
    end
    line1 = sprintf('f := FreeGroup( %s );;', strjoin(cellfun(@(n) ['"' n '"'], names, 'uniform', 0), ', '));
    line2 = sprintf('g := f / [ %s ];;', strjoin(gapRelators, ', '));
    line3 = 'iso := IsomorphismPermGroup(g);;';
    lineRest = arrayfun(@(i) sprintf('ListPerm(Inverse(Image(iso, g.%d)));', i), 1:length(names), 'uniform', 0);
    lines = strjoin(horzcat({line1 line2 line3}, lineRest), '\n');
    tfile = tempname();
    fid = fopen(tfile, 'wt');
    fprintf(fid, lines);
    fclose(fid);
    [status, result] = system([replab.globals.gapBinaryPath ' -q <' tfile]);
    delete(tfile);
    gens = cellfun(@strtrim, strsplit(result, '\n'), 'uniform', 0);
    mask = cellfun(@isempty, gens);
    gens = gens(~mask);
    for i = 1:length(gens)
        gen = strtrim(gens{i});
        assert(gen(1) == '[');
        assert(gen(end) == ']');
        gen = str2num(gen);
        gens{i} = gen;
    end
    ds = max(cellfun(@length, gens));
    for i = 1:length(gens)
        gens{i} = [gens{i} length(gens{i})+1:ds];
    end
end
