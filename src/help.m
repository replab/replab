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
        obj = package.lookupMemberName(restNameParts);
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
    sub = codeBase.subPackages(package.nameParts);
    if ~isempty(sub)
        disp('  Subpackages:');
        for i = 1:length(sub)
            disp(sprintf('    %s', sub{i}));
        end
        disp(' ');
    end
    disp('  Members:')
    fn = fieldnames(package.members);
    for i = 1:length(fn)
        name = fn{i};
        member = package.members.(name);
        if ~isempty(member.docLines)
            disp(sprintf('    %s: %s', member.headerStr, member.docLines{1}));
        else
            disp(sprintf('    %s', member.headerStr));
        end
    end
end

function help_function(codeBase, packageNameParts, fun)
    disp(['In package ' strjoin(packageNameParts, '.')]);
    disp(' ')
    disp(fun.declaration);
    disp(' ');
    for i = 1:length(fun.docLines)
        disp([fun.docLines{i} ' ']);
    end
end
