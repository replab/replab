function c = crawl(rootFolder)
% Crawls the RepLAB source repository
%
% Args:
%   rootFolder (charstring): Absolute path of the source directory (usually '$REPLAB_ROOT/src')
%
% Returns:
%   `.CodeBase`: The parsed code base

    % First, we crawl the source code repository and enumerate files.
    % This is quick.
    toExplore = {{}};
    % toExplore represents a stack of subpaths to explore
    %
    % toExplore is a row cell array, each element inside
    % is a cell array of char strings, which represent a
    % sequence of subfolders of pathStr
    packages = struct;
    % packages is a struct whose field names are pacakge names
    nFiles = 0;
    while length(toExplore) > 0
        % the current path explored
        subpath = toExplore{1};
        packageNameParts = cellfun(@(x) x(2:end), subpath, 'uniform', 0);
        toExplore = toExplore(2:end);
        % the path to explore
        path = fullfile(rootFolder, subpath{:});
        children = dir(path);
        filenames = {}; % discovered files
        for i = 1:length(children)
            name = children(i).name;
            if isequal(name, '.') || isequal(name, '..')
                % do nothing
            elseif children(i).isdir
                % folder
                assert(name(1) == '+', 'We only support crawling subpackages');
                newsubpath = horzcat(subpath, {name});
                toExplore{1,end+1} = newsubpath;
            elseif isequal(name(end-1:end), '.m')
                % is not a folder and has a Matlab file extension
                filenames{1, end+1} = fullfile(rootFolder, subpath{:}, name);
                nFiles = nFiles + 1;
            end
        end
        fn = replab.infra.shm.encode(packageNameParts);
        packages.(fn) = filenames;
    end

    % Now we parse source code data.
    packageData = {};
    ind = 1; % current file index
    fns = fieldnames(packages);
    pb = replab.infra.repl.ProgressBar(nFiles);
    for i = 1:length(fns);
        fn = fns{i};
        filenames = packages.(fn);
        packageNameParts = replab.infra.shm.decode(fn);
        ownFunctions = {};
        ownClasses = {};
        for j = 1:length(filenames)
            filename = filenames{j};
            pb.step(ind, filename);
            data = replab.infra.CodeTokens.fromFile(filename).parse;
            switch class(data)
              case 'replab.infra.FunctionLikeData'
                ownFunctions{1,end+1} = data;
              case 'replab.infra.ClassData'
                ownClasses{1,end+1} = data;
              otherwise
                error(sprintf('Unknown type %s', class(data)));
            end
            ind = ind + 1;
        end
        packageData{1,end+1} = replab.infra.PackageData(packageNameParts, ownFunctions, ownClasses);
    end
    pb.finish;
    c = replab.infra.CodeBase(rootFolder, packageData);
end
