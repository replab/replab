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
        switch length(restNameParts)
          case 0
            % we have a package
            help_package(codeBase, package);
          case 1
            % we have a package member
            memberName = restNameParts{1};
            if ~package.hasMember(memberName)
                  err = 'In package %s, we are looking for member %s, and the member %s is unavailable.';
                  error(sprintf(err, package.fullName, memberName, memberName));
            end
            obj = package.member(memberName);
            if isa(obj, 'replab.infra.Function')
                help_function(codeBase, obj);
            elseif isa(obj, 'replab.infra.Class')
                help_class(codeBase, obj);
            else
                error(sprintf('The object %s should not be directly a member of package %s', obj.headerStr, package.fullName));
            end
          case 2
            % we have a class and a class member
            className = restNameParts{1};
            classElementName = restNameParts{2};
            if ~package.hasMember(className)
                err = 'In package %s, we are looking for class %s with member %s, and the class %s is unavailable.';
                error(sprintf(err, package.fullName, className, classElementName, className));
            end
            class = package.member(className);
            if ~isa(class, 'replab.infra.Class')
                err = 'In package %s, we are looking for class %s with member %s, and %s is not a class.';
                error(sprintf(err, package.fullName, className, classElementName, className));
            end
            if ~class.hasInheritedMember(codeBase, classElementName)
                err = 'In package %s, we are looking for class %s with member %s, and the member %s is not found.';
                error(sprintf(err, package.fullName, className, classElementName, classElementName));
            end
            help_classElement(codeBase, class, classElementName);
          otherwise
            err = 'In package %s, we are looking for the path %s which is too long to describe something sensible';
            error(sprintf(err, package.fullName, strjoin(restNameParts, '.')));
        end
    else
        if isempty(replab.Parameters.matlabHelpPath)
            error('The matlab help path was not captured. Please use replab_addpath first.');
        else
            % We call the matlab help function
            currentPath = strrep(pwd, '\', '/');
            
            if ~replab.platformIsOctave
                cd(replab.Parameters.matlabHelpPath);
                message = [];
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
    disp(['Package ' strjoin(package.nameParts)]);
    disp(' ');
    sub = codeBase.subPackagesNames(package.nameParts);
    if ~isempty(sub)
        disp('  Subpackages:');
        for i = 1:length(sub)
            name = sub{i};
            fullName = strjoin(horzcat(package.nameParts, {name}), '.');
            if replab.platformIsOctave
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
    table = cell(0, 2);
    for i = 1:length(fn)
        name = fn{i};
        member = package.members.(name);
        if replab.platformIsOctave
            ref = name;
        else
            ref = sprintf('<a href="matlab: help(''%s'')">%s</a>', member.fullName, name);
        end
        switch member.kind
          case 'class'
            k = 'cls';
          case 'function'
            k = 'fun';
          otherwise
            k = member.kind;
        end
        table{i,1} = sprintf('    %s (%s)', ref, k);
        table{i,2} = '';
        if ~member.doc.isempty
            table{i,2} = [' ' member.doc.firstLine];
        end
    end
    t = replab.infra.align(table, 'll');
    for i = 1:length(t)
        disp(t{i});
    end
end

function help_class(codeBase, class)
    disp(' ');
    str = ['Class ' class.fullName];
    switch class.nParents 
      case 0
      case 1
        str = [str ' with parent '];
      otherwise
        str = [str ' with parents '];
    end
    sep = '';
    for i = 1:class.nParents
        pn = class.parentName(i);
        str = sprintf('%s%s<a href="matlab: help(''%s'')">%s</a>', str, sep, pn, pn);
        sep = ', ';
    end
    disp(str);
    class.doc.dispFilteredLines;
    disp(' ');
end

function help_classElement(codeBase, class, elementName)
    elements = class.findInheritedMember(codeBase, elementName);
    k = elements{1}.kind;
    disp(sprintf('%s.%s (%s)', class.fullName, elementName, k));
    disp('present in:');
    doc = [];
    for i = 1:length(elements)
        el = elements{i};
        s = strjoin(horzcat(el.packageNameParts, {el.className}), '.');
        if el.doc.isempty
            d = '(no doc)';
        else
            d = '(doc)';
            if isempty(doc)
                doc = el;
                d = '(doc, shown below)';
            end
        end
        disp(sprintf('  %s %s', s, d));
    end
    disp(' ');
    if isequal(k, 'method');
        el1 = elements{1};
        disp(['Declaration in ' el1.classFullName ':']);
        disp(el1.declaration);
        disp(' ');
    end
    if ~isempty(doc)
        if isequal(k, 'method')
            disp(['Documentation in ' doc.classFullName ':']);
            disp(doc.declaration);
        end
        doc.doc.dispFilteredLines;
    end
    disp(' ');    
end

function help_function(codeBase, fun)
    disp(' ');
    disp(['In package ' strjoin(fun.packageNameParts, '.') ':']);
    disp(fun.declaration);
    fun.doc.dispFilteredLines;
    disp(' ');
end
