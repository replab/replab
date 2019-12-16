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
    [srcRoot, name, ~] = fileparts(mfilename('fullpath'));
    % RepLAB root folder
    [root, ~] = fileparts(srcRoot);
    % Folder with tests
    testRoot = fullfile(root, 'tests');
    % Subfolder with doctests
    doctestRoot = fullfile(root, 'tests', 'doctest');

    %% Prepare test directory structure
    switch exist(doctestRoot)
      case 7
        disp('Doctest directory exists, removing it');
        if replab.platformIsOctave
            confirm_recursive_rmdir (false, 'local');
        end
        rmdir(doctestRoot, 's');
      case 0
        disp('Doctest directory does not exist yet');
      otherwise
        error('Unknown type')
    end

    % Create subfolder if inexistent
    [success, message, messageid] = mkdir(testRoot, 'doctest');
    
    disp('Crawling code base');
    codeBase = replab.infra.CodeBase.crawl(fullfile(root, 'src'));
    
    disp('Writing tests');
    codeBase.writeDocTests(doctestRoot);
end
