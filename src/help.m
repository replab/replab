function help(varargin)
% Help function
    persistent codeBase
    if (length(varargin) == 1) && (length(varargin{1}) >= 6) && (isequal(varargin{1}(1:6), 'replab'))
        name = varargin{1};
        if isempty(codeBase)
            [srcRoot, ~, ~] = fileparts(mfilename('fullpath'));
            codeBase = replab.infra.CodeBase.crawl(srcRoot);
        end
        parts = strsplit(name, '.');
        [package packageNameParts restNameParts] = codeBase.lookupGreedy(parts);
        try
            obj = package.lookupMemberName(restNameParts);
        catch
        end
        if isa(obj, 'replab.infra.Function')
            help_function(codeBase, packageNameParts, obj);
        elseif isa(obj, 'replab.infra.Class')
            help_class(codeBase, packageNameParts, obj);
        elseif isa(obj, 'replab.infra.Package')
            help_package(codeBase, obj);
        end
    else
        if isempty(replab.Parameters.matlabHelpPath)
            error('The matlab help path was not captured. Please use replab_addpath first.');
        else
            % We call the matlab help function
            currentPath = strrep(pwd, '\', '/');
            
            isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
            if ~isOctave
                cd(replab.Parameters.matlabHelpPath);
                try
                    help(varargin{:});
                catch message
                end
                cd(currentPath);
                if ~isempty(message)
                    error(message);
                end
            else
                % In some versions of octave earlier than 5.1.0, the
                % current path had a lower priority than the path order.
                % Then we also need replab's path...
                
                replabHelpPath = fileparts(which('replab_addpaths'));
                replabHelpPath = [strrep(replabHelpPath, '\', '/'), '/src'];
                
                cd(replab.Parameters.matlabHelpPath);
                addpath(replab.Parameters.matlabHelpPath);
                message = [];
                try
                    help(varargin{:});
                catch message
                end
                cd(currentPath);
                addpath(replabHelpPath);
                if ~isempty(message)
                    error(message);
                end
            end
        end
    end
end

function help_package(codeBase, package)
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

    disp(['Package ' strjoin(package.nameParts)]);
    disp(' ');
    sub = codeBase.subPackagesNames(package.nameParts);
    if ~isempty(sub)
        disp('  Subpackages:');
        for i = 1:length(sub)
            name = sub{i};
            fullName = strjoin(horzcat(package.nameParts, {name}), '.');
            if isOctave
                ref = name;
            else
                ref = sprintf('<a href="matlab: help(''%s'')">%s</a>', fullName, name);
            end
            disp(sprintf('    %s', ref));
        end
        disp(' ');
    end
    disp('  Members:')
    fn = fieldnames(package.members);
    table = cell(0, 3);
    for i = 1:length(fn)
        name = fn{i};
        member = package.members.(name);
        if isOctave
            ref = name;
        else
            ref = sprintf('<a href="matlab: help(''%s'')">%s</a>', member.fullName, name);
        end
        table{i,1} = ['    ' ref];
        table{i,2} = [' ' member.kind];
        table{i,3} = '';
        if ~member.doc.isempty
            table{i,3} = [' ' member.doc.firstLine];
        end
    end
    t = replab.infra.align(table, 'lll');
    for i = 1:length(t)
        disp(t{i});
    end
end

function help_class(codeBase, class)
end

function help_method(codeBase, method)
end

function help_property(codeBase, property)
end

function help_function(codeBase, fun)
    disp(['In package ' strjoin(packageNameParts, '.')]);
    disp(' ')
    disp(fun.declaration);
    disp(' ');
    for i = 1:length(fun.docLines)
        disp([fun.docLines{i} ' ']);
    end
end
