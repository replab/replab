function help(varargin)
% Help function overloading

    % extract own name
    [pathStr, helpFunctionName, extension] = fileparts(which(mfilename));

    persistent codeBase
    
    if (length(varargin) == 1) && ischar(varargin{1})
        % we check if the the subject starts with a variable
        dotPositions = strfind(varargin{1}, '.');
        if isempty(dotPositions)
            firstPart = varargin{1};
            secondPart = '';
        else
            firstPart = varargin{1}(1:dotPositions(1)-1);
            secondPart = varargin{1}(dotPositions(1)+1:end);
        end
        isObject = false;
        try
            isObject = evalin('caller', ['isobject(', firstPart, ')']);
        catch
        end
        if isObject
            % If asking help regarding an object, we see if it is a replab
            % object
            type = evalin('caller', ['class(', firstPart, ')']);
            if (length(type) >= 6) && isequal(type(1:6), 'replab')
                if isempty(secondPart)
                    varargin{1} = type;
                else
                    varargin{1} = [type, '.', secondPart];
                end
            end
        end
    end
    
    if (length(varargin) == 1) && (length(varargin{1}) >= 6) && (isequal(varargin{1}(1:6), 'replab'))
        % We are looking for a replab-related help
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
            help_package(codeBase, package, helpFunctionName);
          case 1
            % we have a package member
            memberName = restNameParts{1};
            if ~package.hasMember(memberName)
                  err = 'In package %s, we are looking for member %s, and the member %s is unavailable.';
                  error(sprintf(err, package.fullName, memberName, memberName));
            end
            obj = package.member(memberName);
            if isa(obj, 'replab.infra.Function')
                help_function(codeBase, obj, helpFunctionName);
            elseif isa(obj, 'replab.infra.Class')
                help_class(codeBase, obj, helpFunctionName);
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
            classe = package.member(className);
            if ~isa(classe, 'replab.infra.Class')
                err = 'In package %s, we are looking for class %s with member %s, and %s is not a class.';
                error(sprintf(err, package.fullName, className, classElementName, className));
            end
            if ~classe.hasInheritedMember(codeBase, classElementName)
                err = 'In package %s, we are looking for class %s with member %s, and the member %s is not found.';
                error(sprintf(err, package.fullName, className, classElementName, classElementName));
            end
            help_classElement(codeBase, classe, classElementName, helpFunctionName);
          otherwise
            err = 'In package %s, we are looking for the path %s which is too long to describe something sensible';
            error(sprintf(err, package.fullName, strjoin(restNameParts, '.')));
        end
    else
        if isempty(replab.Parameters.matlabHelpPath)
            error('The matlab help path was not captured. Please use replab_addpaths first.');
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

function help_package(codeBase, package, helpFunctionName)
    keyword = strjoin(package.nameParts);
    replab.infra.dispH(['Package ' strjoin(package.nameParts), ':'], keyword, helpFunctionName);
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
                ref = sprintf('<a href="matlab: %s(''%s'')">%s</a>', helpFunctionName, fullName, name);
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
            ref = sprintf('<a href="matlab: %s(''%s'')">%s</a>', helpFunctionName, member.fullName, name);
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

function help_class(codeBase, class, helpFunctionName)
    keyword = class.name;
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
        str = sprintf('%s%s%s', str, sep, pn);
        sep = ', ';
    end
    replab.infra.dispH(str, keyword, helpFunctionName);
    class.doc.dispFilteredLines(keyword, helpFunctionName);
    disp(' ');
end

function help_classElement(codeBase, class, elementName, helpFunctionName)
    elements = class.findInheritedMember(codeBase, elementName);
    k = elements{1}.kind;
    keyword = elementName;
    replab.infra.dispH(sprintf('%s.%s (%s)', class.fullName, elementName, k), keyword, helpFunctionName);
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
        replab.infra.dispH(sprintf('  %s %s', s, d), keyword, helpFunctionName);
    end
    disp(' ');
    if isequal(k, 'method')
        el1 = elements{1};
        replab.infra.dispH(['Declaration in ' el1.classFullName ':'], keyword, helpFunctionName);
        replab.infra.dispH(el1.declaration, keyword, helpFunctionName);
        disp(' ');
    end
    if ~isempty(doc)
        if isequal(k, 'method')
            replab.infra.dispH(['Documentation in ' doc.classFullName ':'], keyword, helpFunctionName);
            replab.infra.dispH(doc.declaration, keyword, helpFunctionName);
        end
        doc.doc.dispFilteredLines(keyword, helpFunctionName);
    end
    disp(' ');    
end

function help_function(codeBase, fun, helpFunctionName)
    keyword = fun.name;
    disp(' ');
    replab.infra.dispH(['In package ' strjoin(fun.packageNameParts, '.') ':'], keyword, helpFunctionName);
    disp(' ');
    replab.infra.dispH(fun.declaration, keyword, helpFunctionName);
    fun.doc.dispFilteredLines(keyword, helpFunctionName);
    disp(' ');
    if ~replab.platformIsOctave
        disp(['    <a href="matlab: rdoc(''', strjoin(fun.packageNameParts, '.'), ''')">Reference page for ', strjoin(fun.packageNameParts, '.'), '</a>'])
    end
end
