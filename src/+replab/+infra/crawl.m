function c = crawl(rootFolder)
% Crawls the RepLAB source repository
%
% Args:
%   rootFolder (charstring): Absolute path of the source directory (usually '$REPLAB_ROOT/src')
%
% Returns:
%   `.CodeBase`: The parsed code base
    packageData = {};
    % toExplore represents a stack of subpaths to explore
    % toExplore is a row cell array, each element inside
    % is a cell array of char strings, which represent a
    % sequence of subfolders of pathStr
    toExplore = {{}};
    while length(toExplore) > 0
        % the current path explored
        subpath = toExplore{1};
        packageNameParts = cellfun(@(x) x(2:end), subpath, 'uniform', 0);
        toExplore = toExplore(2:end);
        % the path to explore
        path = fullfile(rootFolder, subpath{:});
        children = dir(path);
        ownFunctions = {};
        ownClasses = {};
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
                filename = fullfile(rootFolder, subpath{:}, name);
                data = replab.infra.CodeTokens.fromFile(filename).parse;
                switch class(data)
                  case 'replab.infra.FunctionLikeData'
                    ownFunctions{1,end+1} = data;
                  case 'replab.infra.ClassData'
                    ownClasses{1,end+1} = data;
                  otherwise
                    error(sprintf('Unknown type %s', class(data)));
                end
            end
        end
        packageData{1,end+1} = replab.infra.PackageData(packageNameParts, ownFunctions, ownClasses);
    end
    c = replab.infra.CodeBase(rootFolder, packageData);
end
