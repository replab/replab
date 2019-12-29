function writeDocTests(doctestPath, el)
% Writes the doctests for the given package element
%
% Args:
%   doctestPath (charstring): Base folder for doctest generation, not including trailing path separator
%                             That folder must exist.
%   el (`replab.infra.SourceElement`): Element to inspect for doctests
    doctests = replab.infra.DocTest.parseDoc(el.doc);
    lineNumbers = cellfun(@(x) 2, doctests);
    elementNames = cellfun(@(x) '', doctests, 'uniform', 0);
    testNumbers = 1:length(doctests);
    if isa(packageElement, 'replab.infra.Class')
        memberNames = packageElement.memberNames;
        for i = 1:length(memberNames)
            member = packageElement.member(memberNames{i});
            newDoctests = replab.infra.DocTest.parseDoc(member.doc);
            for j = 1:length(newDoctests)
                doctests{1,end+1} = newDoctests{j};
                lineNumbers(1,end+1) = member.lineNumber + 1;
                elementNames{1,end+1} = memberNames{i};
                testNumbers(1,end+1) = j;
            end
        end
    end
    if ~isempty(doctests)
        base = doctestPath;
        for i = 1:length(packageElement.packageNameParts)
            part = packageElement.packageNameParts{i};
            switch exist(fullfile(base, part))
              case 0
                  assert(mkdir(base, part), sprintf('Could not create %s in directory %s', part, base));
              case 7
                % directory already exists
              otherwise
                error(sprintf('%s already exists in directory %s, but is not a folder', part, base));
            end
            base = fullfile(base, part);
        end
        functionName = [packageElement.name 'Test'];
        filename = fullfile(base, [functionName '.m']);
        fid = fopen(filename, 'w');
        fprintf(fid, 'function test_suite = %s()\n', functionName);
        fprintf(fid, '  disp([''Setting up tests in '', mfilename()]);\n');
        fprintf(fid, '  try\n');
        fprintf(fid, '    test_functions = localfunctions();\n');
        fprintf(fid, '  catch\n');
        fprintf(fid, '  end\n');
        fprintf(fid, '  initTestSuite;\n');
        fprintf(fid, 'end\n\n');
        for i = 1:length(doctests)
            replab.infra.writeDocTest(fid, packageElement.fullFilename, lineNumbers(i), elementNames{i}, testNumbers(i), doctests{i});
        end
        if length(doctests) == 1
            disp('  One doctest found.');
        else
            fprintf('  %d doctests found.\n', length(doctests));
        end
        fclose(fid);
    end
end
