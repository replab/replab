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
    % Create subfolder if inexistent
    [success, message, messageid] = mkdir(testRoot, 'doctest');
    
    %% Recurse source directory
    codeFiles = {};
    % toExplore represents a stack of subpaths to explore
    % toExplore is a row cell array, each element inside
    % is a cell array of char strings, which represent a
    % sequence of subfolders of pathStr
    toExplore = {{}};
    subPaths = {};
    while length(toExplore) > 0
        subpath = toExplore{1};
        toExplore = toExplore(2:end);
        % the path to explore
        path = fullfile(srcRoot, subpath{:});
        children = dir(path);
        for i = 1:length(children)
            name = children(i).name;
            if isequal(name, '.') || isequal(name, '..')
                % do nothing
            elseif children(i).isdir
                % folder
                newsubpath = horzcat(subpath, name);
                toExplore{end+1} = newsubpath;
                subPaths{end+1} = newsubpath;
            elseif isequal(name(end-1:end), '.m')
                % is not a folder and has a Matlab file extension
                codeFiles{end+1} = horzcat(subpath, name);
            end
        end
    end
    disp(sprintf('Found %d code files in %d folders', length(codeFiles), length(subPaths)));
    
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
    disp('Creating subfolders in the test directory');
    for i = 1:length(subPaths)
        subpath = subPaths{i};
        testsubpath = cellfun(@(s) strrep(s, '+', ''), subpath, 'uniform', 0);
        parent = fullfile(doctestRoot, testsubpath{1:end-1});
        new = testsubpath{end};
        [success, message, messageid] = mkdir(parent, new);
    end
    
    %% Parse source code files looking for doctest blocks
    for i = 1:length(codeFiles)
        subpath = codeFiles{i};
        testsubpath = cellfun(@(s) strrep(s, '+', ''), subpath, 'uniform', 0);
        sourceFile = fullfile(srcRoot, subpath{:});
        fprintf('Processing %s', sourceFile);
        testsubpath{end} = strrep(testsubpath{end}, '.m', 'Test.m');
        testFile = fullfile(doctestRoot, testsubpath{:});
        fid = fopen(sourceFile, 'r');
        l = fgetl(fid);
        % we read the file line by line, with a state machine
        % state = 0: waiting for a >>> prompt (corresponding to a new test)
        % state = 1: reading expected output lines
        % state = 2: reading expected output, already got one blank line
        tests = {};
        test = {};
        subtestNb = 1;
        messages = {};
        message = {};
        state = 0;
        lineCount = 1;
        while ~isequal(l, -1)
            l = strtrim(l);
            if length(l) < 1 || l(1) ~= '%'
                % no longer a comment block
                if numel(test) > 1
                    % if there is a parsed test, register it
                    tests{end+1} = test;
                    test = {};
                    messages{end+1} = message;
                    message = {};
                    subtestNb = 1;
                end
                state = 0;
            else
                % we have a comment line, remove the leading % and trim
                l = strtrim(l(2:end));
                if length(l) > 3 && isequal(l(1:3), '>>>')
                    % register current test if there was one and new
                    % command belongs to a new test
                    if (numel(test) >= 1)
                        if state == 0
                            % command belongs to a new test
                            tests{end+1} = test;
                            test = {};
                            messages{end+1} = message;
                            message = {};
                            subtestNb = 1;
                        else
                            % command adds to current test
                            subtestNb = subtestNb + 1;
                        end
                    end
                    test{subtestNb} = {strtrim(l(4:end))};
                    message{subtestNb} = ['Doctest failure on line ', num2str(lineCount), ' of ', sourceFile];
                    state = 1;
                elseif length(l) == 0
                    switch state
                      case 0 % nothing
                      case 1
                        state = 2;
                      case 2
                        if numel(test) > 1
                            tests{end+1} = test;
                            messages{end+1} = message;
                        end
                        test = {};
                        message = {};
                        subtestNb = 1;
                        state = 0;
                    end
                else
                    if state > 0
                        test{subtestNb}{1, end+1} = l;
                    end
                end
            end
            l = fgetl(fid);
            lineCount = lineCount + 1;
        end
        fclose(fid);
        fprintf('-- Read %d tests\n', length(tests));
        if length(tests) > 0
            % If there are tests in the source file, write them
            fid = fopen(testFile, 'w');
            functionName = strrep(testsubpath{end}, '.m', '');
            fprintf(fid, 'function test_suite = %s()\n', functionName);
            fprintf(fid, '  disp([''Setting up tests in '', mfilename()]);\n');
            fprintf(fid, '  try\n');
            fprintf(fid, '    test_functions = localfunctions();\n');
            fprintf(fid, '  catch\n');
            fprintf(fid, '  end\n');
            fprintf(fid, '  initTestSuite;\n');
            fprintf(fid, 'end\n\n');
            for i = 1:length(tests)
                test = tests{i};
                message = messages{i};
                fprintf(fid, 'function test%d\n', i);
                for j = 1:length(test)
                    outs = cellfun(@(x) ['''' strrep(x, '''', '''''') ''''], test{j}(2:end), 'uniform', 0);
                    out = strjoin(outs, ' ');
                    fprintf(fid, '  out = evalc(''%s'');\n', strrep(test{j}{1}, '''', ''''''));
                    fprintf(fid, '  assertEqualEvalcOutput(out, {%s}, ''%s'');\n', out, message{j});
                end
                fprintf(fid, 'end\n\n');
            end
            fclose(fid);
        end
    end
end
