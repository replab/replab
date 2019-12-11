function help(name)
% Help function
    persistent codeBase
    if isequal(name(1:6), 'replab')
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
        error('JDB: put your special code here');
    end
end

function help_package(codeBase, package)
    disp(['Package ' strjoin(package.nameParts)]);
    disp(' ');
    sub = codeBase.subPackagesNames(package.nameParts);
    if ~isempty(sub)
        disp('  Subpackages:');
        for i = 1:length(sub)
            name = sub{i};
            fullName = strjoin(horzcat(package.nameParts, {name}), '.');
            ref = sprintf('<a href="matlab: help(''%s'')">%s</a>', fullName, name);
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
        ref = sprintf('<a href="matlab: help(''%s'')">%s</a>', member.fullName, name);
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
