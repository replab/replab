function writeElementDocTests(doctestPath, el)
% Write the doctests for the given source element (class or function)
%
% Args:
%   doctestPath (charstring): Base folder for doctest generation, not including trailing path separator
%                             That folder must exist.
%   el (`+replab.+infra.SourceElement`): Element to inspect for doctests
    s = replab.infra.doctests.namedTestsInElement(el);
    if isempty(fieldnames(s))
        return
    end
    folder = fullfile(doctestPath, el.package.packagePath{:});
    % Make path
    switch exist(folder)
      case 0
        assert(mkdir(folder), sprintf('Could not create directory %s', folder));
      case 7
        % directory already exists
      otherwise
        error(sprintf('%s already exists but is not a folder', folder));
    end
    testFunName = [el.name 'Test'];
    testFilename = fullfile(folder, [testFunName '.m']);
    fid = fopen(testFilename, 'w');
    initFN = fullfile(el.codeBase.rootFolder, '+replab', '+infra', '+doctests', 'doctest_init.liquid');
    init = replab.lobster.Template.load(initFN);
    filename = el.absoluteFilename;
    fprintf(fid, init.render(struct('testFunName', testFunName)));
    names = fieldnames(s);
    for i = 1:length(names)
        name = names{i};
        values = s.(name);
        if length(values) == 1
            innerNames = {name};
        else
            innerNames = arrayfun(@(j) sprintf('%s%d', name, j), 1:length(values), 'uniform', 0);
        end
        for j = 1:length(values)
            fprintf(fid, 'function test_%s\n', innerNames{j});
            fprintf(fid, '  filename = ''%s'';\n', filename);
            v = values{j};
            for k = 1:v.nCommands
                fprintf(fid, '  out = evalc(%s);\n', v.quotedCommand(k));
                fprintf(fid, '  assertEqualEvalcOutput(out, %s, filename, %d);\n', v.quotedOutput(k), v.lineNumbers(k));
            end
            fprintf(fid, 'end\n');
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
end
