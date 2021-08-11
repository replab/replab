function writeNotebook(notebooksPath, notebook)
% Write the notebook code into
%
% Args:
%   doctestPath (charstring): Base folder to hold the notebook, not including trailing path separator
%                             That folder must exist.
%   notebook (cell array of charstring): notebook path in three parts

    folder = fullfile(notebooksPath, notebook{2});
    % Make path
    switch exist(folder)
      case 0
        assert(mkdir(folder), sprintf('Could not create directory %s', folder));
      case 7
        % directory already exists
      otherwise
        error(sprintf('%s already exists but is not a folder', folder));
    end

    testFunName = [notebook{3}(1:end-2) 'Test'];
    testFilename = fullfile(folder, [testFunName '.m']);
    filenameIn = fullfile(notebook{1}, notebook{2}, notebook{3});
    fidIn = fopen(filenameIn, 'r');
    fidOut = fopen(testFilename, 'w');
    initFN = fullfile(replab.globals.replabPath, 'src', '+replab', '+infra', '+doctests', 'doctest_init.liquid');
    init = replab.lobster.Template.load(initFN);
    fprintf(fidOut, init.render(struct('testFunName', testFunName)));
    
    fprintf(fidOut, 'function test_%s\n', testFunName);
    fprintf(fidOut, '  filename = ''%s'';\n', filenameIn);
    % Copy file over
    line = fgetl(fidIn);
    while ischar(line)
        % We copy the line over
    	if (length(line) >= 3) && isequal(line(end-2:end), '...')
            fprintf(fidOut, '%s\n', line);
    	else
            % If the line finishes, we add an extra ';' to limit the verbose level
            fprintf(fidOut, '%s;\n', line);
    	end
        line = fgetl(fidIn);
    end
    fprintf(fidOut, 'end\n');
    fprintf(fidOut, '\n');

    fclose(fidOut);
end
