function replab_generatedoctests
% Extracts doctests from the source files in the src directory
%
% The code recurses the source files, and looks for comment lines starting with ``>>>``.
% The text following ``>>>`` is interpreted as a command. The lines following ``>>>`` are
% interpreted as the test expected output, when run in the command line. Expected test output
% ends 1) when the next line is no longer a comment line, 2) when another comment line starts
% with ``>>>``, 3) when two blank comment lines are provided.
%
% In expected test output, blank lines are ignored. Errors are not handled.
%
% Test files are written in the MoXUnit format, using functions and subfunctions as the
% rest of the RepLAB test suite; the helper function ``assertEqualEvalcOutput`` is used
% to verify the test output.
    rp = replab.settings.replabPath;
    srcRoot = fullfile(rp, 'src');
    testRoot = fullfile(rp, 'tests');
    doctestRoot = fullfile(rp, 'tests', 'doctest');

    %% Prepare test directory structure
    switch exist(doctestRoot)
      case 7
        disp('Doctest directory exists, removing it');
        replab.compat.rmdirRec(doctestRoot);
      case 0
        disp('Doctest directory does not exist yet');
      otherwise
        error('Unknown type')
    end

    % Create subfolder if inexistent
    [success, message, messageid] = mkdir(testRoot, 'doctest');

    disp('Crawling code base');
    cb = replab.infra.crawl(srcRoot);
    disp('Generating doctests');
    af = cb.allFunctions;
    for i = 1:length(af)
        replab.infra.doctests.writeElementDocTests(doctestRoot, af{i});
    end
    ac = cb.allClasses;
    for i = 1:length(ac)
        replab.infra.doctests.writeElementDocTests(doctestRoot, ac{i});
    end
end
