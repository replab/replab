function codeFiles = replab_generate_doctests
    [srcRoot, name, ~] = fileparts(mfilename('fullpath'));
    % RepLAB root folder
    [root, ~] = fileparts(srcRoot);
    % Folder with tests
    testRoot = fullfile(root, 'tests');
    % Subfolder with doctests
    doctestRoot = fullfile(root, 'tests', 'doctest');
    % Create subfolder if inexistent
    [success, message, messageid] = mkdir(testRoot, doctestRoot);
    % Recurse source directory
    codeFiles = {};
    % subpaths represents a stack of subpaths to explore
    % subpaths is a row cell array, each element inside
    % is a cell array of char strings, which represent a
    % sequence of subfolders of pathStr
    % each element in 
    subpaths = {{}};
    while length(subpaths) > 0
        subpath = subpaths{end};
        subpaths = subpaths(1:end-1);
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
                subpaths{end+1} = newsubpath;
            elseif isequal(name(end-1:end), '.m')
                % is not a folder and has a Matlab file extension
                codeFiles{end+1} = horzcat(subpath, name);
            end
        end
    end
    for i = 1:length(codeFiles)
        subpath = codeFiles{i};
        testsubpath = cellfun(@(s) strrep(s, '+', ''), subpath, 'uniform', 0);
        for j = 1:length(subpath)-1
            parent = fullfile(doctestRoot, testsubpath{1:j-1});
            new = testsubpath{j};
            [success, message, messageid] = mkdir(parent, new);
        end
    end
end

