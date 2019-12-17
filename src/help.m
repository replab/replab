function help(varargin)
% Displays help on a function/object/...
%
% Arg:
%     varargin: - if a single element, provides help for this element,
%                 whether it is a variable, object, function, class, 
%                 package, ...
%               - if two elements, with one of them being the string '-f'
%                 or '--full', provides full help on this element
%
% Note:
%     This function overloads matlab's 'help' function. Renaming this file
%     does not alter its functionality.

    % extract own name
    [pathStr, helpFunctionName, extension] = fileparts(which(mfilename));

    persistent codeBase
    
    % Are we in full mode or not?
    flags = {};
    if (length(varargin) == 2)
        if isequal(varargin{1}, '-f') || isequal(varargin{1}, '--full')
            flags = {'-f'};
            varargin = {varargin{2}};
        elseif isequal(varargin{2}, '-f') || isequal(varargin{2}, '--full')
            flags = {'-f'};
            varargin = {varargin{1}};
        end
    end
    
    % Are we asking help on an object or method applied to an object?
    if (length(varargin) == 1) && ischar(varargin{1})
        % we check if the subject starts with a variable
        dotPositions = strfind(varargin{1}, '.');
        if isempty(dotPositions)
            firstPart = varargin{1};
            secondPart = '';
        else
            firstPart = varargin{1}(1:dotPositions(1)-1);
            secondPart = varargin{1}(dotPositions(1)+1:end);
        end
        isVar = 0;
        try
            isVar = evalin('caller', ['exist(''', firstPart, ''')']);
        catch
        end
        if isVar == 1
            % object of interest is a variable or object (or method applied
            % to an object)
            type = evalin('caller', ['class(', firstPart, ')']);
            if (length(type) >= 6) && isequal(type(1:6), 'replab')
                % object is a replab object
                if isempty(secondPart)
                    varargin{1} = type;
                    replab.infra.dispH([firstPart, ' is an object of type ', type, '.'], firstPart, helpFunctionName, fullMode);
                    disp(' ');
                else
                    varargin{1} = [type, '.', secondPart];
                end
                replab.infra.dispH(['--- help for ', varargin{1}, ' ---'], varargin{1}, helpFunctionName, fullMode);
            else
                % non-replab variable/object. We need to make a copy to our
                % workspace for the next function to see it
                eval([firstPart, ' = evalin(''caller'', ''', firstPart, ''');']);
            end
        end
    end
    
    if (length(varargin) == 1) && (length(varargin{1}) >= 6) && (isequal(varargin{1}(1:6), 'replab'))
        % We are looking for a replab-related help
        name = varargin{1};
        if isempty(codeBase)
            fprintf('Building help index...');
            [srcRoot, ~, ~] = fileparts(mfilename('fullpath'));
            codeBase = replab.infra.CodeBase.crawl(srcRoot);
            disp('done.');
        end
        parts = strsplit(name, '.');
        el = codeBase.get(parts{:});
        switch class(el)
          case 'replab.infra.Package'
            help_package(el, helpFunctionName, flags);
          case 'replab.infra.Class'
            help_class(codeBase, el, helpFunctionName, flags);
          case 'replab.infra.InheritedClassElement'
            % help_classElement(codeBase, classe, classElementName, helpFunctionName, fullMode);
          case 'replab.infra.ConcreteClassElement'
            % help_classElement(codeBase, classe, classElementName, helpFunctionName, fullMode);
          case 'replab.infra.Function'
            help_function(codeBase, el, helpFunctionName, flags);
          otherwise
            error(sprintf('Unknown code base element type %s', class(el)));
        end
    else
        if isempty(replab.Parameters.matlabHelpPath)
            error('The matlab help path was not captured. Please use replab_init first.');
        else
            % We call matlab's help function
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
                
                replabHelpPath = fileparts(which('replab_init'));
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

function help_package(package, helpFunctionName, flags)
    fullName = package.fullIdentifier;
    replab.infra.dispH([' Package ' fullName, ':'], fullName, helpFunctionName, fullMode);
    disp(' ');
    
    % calling convention linkHelp(helpFunctionName, 'Group', 'replab.Group', '-f')
    
    % Package:
    parent = package.parent;
    if ~isempty(parent) && ~isempty(parent.name)
        disp('   In package:');
        packageName = parent.fullIdentifier;
        if replab.platformIsOctave  || (~usejava('desktop'))
            disp(['     ', packageName]);
        else
            dispH(['     ' replab.infra.linkHelp(helpFunctionName, packageName, packageName, flags)]);
        end
        disp(' ');
    end
    
    % Subpackages:
    sub = package.ownSubpackages;
    if ~isempty(sub)
        disp('   Subpackages:');
        for i = 1:length(sub)
            subi = sub{i};
            fullName = subi.fullIdentifier;
            if replab.platformIsOctave || (~usejava('desktop'))
                ref = subi.name;
            else
                if fullMode
                    ref = sprintf('<a href="matlab: %s(''-f'',''%s'')">%s</a>', helpFunctionName, fullName, subi.name);
                else
                    ref = sprintf('<a href="matlab: %s(''%s'')">%s</a>', helpFunctionName, fullName, subi.name);
                end
            end
            disp(sprintf('     %s', ref));
        end
        disp(' ');
    end
    
    % Members:
    disp('   Members:')
    fn = fieldnames(package.nspElements);
    table = cell(length(fn), 2);
    for i = 1:length(fn)
        name = fn{i};
        element = package.nspElements.(name);
        if replab.platformIsOctave || (~usejava('desktop'))
            ref = name;
        else
            if fullMode
                ref = sprintf('<a href="matlab: %s(''-f'',''%s'')">%s</a>', helpFunctionName, element.fullIdentifier, name);
            else
                ref = sprintf('<a href="matlab: %s(''%s'')">%s</a>', helpFunctionName, element.fullIdentifier, name);
            end
        end
        switch class(element)
          case 'replab.infra.Class'
            k = 'cls';
          case 'replab.infra.Function'
            k = 'fun';
          otherwise
            error(sprintf('Unknown type %s', class(element)));
        end
        table{i,1} = sprintf('     %s (%s)', ref, k);
        table{i,2} = '';
        if ~element.doc.isempty
            table{i,2} = [' ' element.doc.firstLine];
        end
    end
    t = replab.infra.align(table, 'll');
    for i = 1:length(t)
        disp(t{i});
    end
    disp(' ');
end


function help_class(codeBase, class, helpFunctionName, fullMode)
    str = [' Class ' class.fullName];
    fullName = class.fullName;

    if ~fullMode
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
        replab.infra.dispH(str, fullName, helpFunctionName, fullMode);
        
        % Doc
        class.doc.dispFilteredLines(fullName, helpFunctionName, fullMode);
        disp(' ');
        
        % Link to ref
        if ~(replab.platformIsOctave  || (~usejava('desktop')))
            disp(['   <a href="matlab: ', helpFunctionName, '(''-f'',''', class.fullName, ''')">Reference page for ', class.fullName, '</a>'])
        end
        disp(' ');
    else
        replab.infra.dispH([str, ':'], fullName, helpFunctionName, fullMode);
        disp(' ');
        
        % Package:
        if ~isempty(class.packageNameParts)
            disp('   In package:');
            packageName = strjoin(class.packageNameParts, '.');
            if replab.platformIsOctave  || (~usejava('desktop'))
                disp(['     ', packageName]);
            else
                disp(['     <a href="matlab: ', helpFunctionName, '(''-f'',''', packageName, ''')">', packageName, '</a>']);
            end
            disp(' ');
        end
        
        % Parent classes:
        if class.nParents > 0
            switch class.nParents
              case 1
                disp('   Parent:');
              otherwise
                disp('   Parents:');
            end
            table = cell(0,2);
            for i = 1:class.nParents
                table{end+1,1} = ['     ', class.parentName(i), ' '];
                [package packageNameParts restNameParts] = codeBase.lookupGreedy(class.parentNameParts(i));
                if package.hasMember(restNameParts{1})
                    table{end,2} = package.member(restNameParts{1}).doc.firstLine;
                else
                    table{end,2} = '';
                end
            end
            t = replab.infra.align(table, 'll');
            for i = 1:length(t)
                replab.infra.dispH(t{i}, fullName, helpFunctionName, fullMode);
            end
            disp(' ');
        end

        % Children classes:
        children = class.childrenNames(codeBase);
        if ~isempty(children)
            switch length(children)
                case 1
                    disp('   Child:');
                otherwise
                    disp('   Children:');
            end
            table = cell(0,2);
            for i = 1:length(children)
                table{end+1,1} = ['     ', children{i}, ' '];
                [package packageNameParts restNameParts] = codeBase.lookupGreedy(strsplit(children{i},'.'));
                if package.hasMember(restNameParts{1})
                    table{end,2} = package.member(restNameParts{1}).doc.firstLine;
                else
                    table{end,2} = '';
                end
            end
            t = replab.infra.align(table, 'll');
            for i = 1:length(t)
                replab.infra.dispH(t{i}, fullName, helpFunctionName, fullMode);
            end
            disp(' ');
        end
        
        % Properties and Methods:
        fn = fieldnames(class.members);
        table = cell(length(fn), 2);
        tableProperties = cell(0,2);
        tableMethods = cell(0,2);
        tableElse = cell(0,2);
        for i = 1:length(fn)
            name = fn{i};
            member = class.members.(name);
            if replab.platformIsOctave || (~usejava('desktop'))
                ref = name;
            else
                ref = sprintf('<a href="matlab: %s(''-f'',''%s'')">%s</a>', helpFunctionName, member.fullName, name);
            end
            switch member.kind
              case 'class'
                tableElse{end+1,1} = sprintf('     %s (%s)', ref, 'cls');
                tableElse{end,2} = '';
                if ~member.doc.isempty
                    tableElse{end,2} = [' ' member.doc.firstLine];
                end
              case 'function'
                tableElse{end+1,1} = sprintf('     %s (%s)', ref, 'fun');
                tableElse{end,2} = '';
                if ~member.doc.isempty
                    tableElse{end,2} = [' ' member.doc.firstLine];
                end
              case 'property'
                tableProperties{end+1,1} = sprintf('     %s', ref);
                tableProperties{end,2} = '';
                if ~member.doc.isempty
                    tableProperties{end,2} = [' ' member.doc.firstLine];
                end
              case 'method'
                tableMethods{end+1,1} = sprintf('     %s', ref);
                tableMethods{end,2} = '';
                if ~member.doc.isempty
                    tableMethods{end,2} = [' ' member.doc.firstLine];
                end
              otherwise
                tableElse{end+1,1} = sprintf('     %s (%s)', ref, member.kind);
                tableElse{end,2} = '';
                if ~member.doc.isempty
                    tableElse{end,2} = [' ' member.doc.firstLine];
                end
            end
        end
        if ~isempty(tableProperties)
            disp('   Properties:');
            t = replab.infra.align(tableProperties, 'll');
            for i = 1:length(t)
                disp(t{i});
            end
            disp(' ');
        end
        if ~isempty(tableMethods)
            disp('   Methods:')
            t = replab.infra.align(tableMethods, 'll');
            for i = 1:length(t)
                disp(t{i});
            end
            disp(' ');
        end
        if ~isempty(tableElse)
            disp('   Additional members:')
            t = replab.infra.align(tableElse, 'll');
            for i = 1:length(t)
                disp(t{i});
            end
            disp(' ');
        end
        
        % Link to help
        if ~(replab.platformIsOctave  || (~usejava('desktop')))
            disp(['   <a href="matlab: ', helpFunctionName, '(''', class.fullName, ''')">Help page for ', class.fullName, '</a>'])
        end
        disp(' ');
    end
    
    % Link to source code
    disp(sprintf('   <a href="matlab: opentoline(''%s'',%d,1)">See source</a>', class.fullFilename, 1));
    disp(' ');
end

function help_classElement(codeBase, class, elementName, helpFunctionName, fullMode)
    elements = class.findInheritedMember(codeBase, elementName);
    k = elements{1}.kind;
    fullName = [class.fullName, '.', elementName];

    % We look for the best declaration and doc
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
    end
    
    if ~fullMode
        replab.infra.dispH([' ', upper(k(1)), k(2:end), ' ', fullName], fullName, helpFunctionName, fullMode);
        
        % Declaration
        if isequal(k, 'method')
            el1 = elements{1};
            disp(' ');
            replab.infra.dispH(['   ', el1.declaration], fullName, helpFunctionName, fullMode);
            disp(' ');
        end
        
        % Doc
        if ~isempty(doc)
            if ~isequal(doc.classFullName, class.fullName)
                replab.infra.dispH(['  Documentation from ', doc.classFullName, ':'], fullName, helpFunctionName, fullMode);
            end
            doc.doc.dispFilteredLines(fullName, helpFunctionName, fullMode);
        end
        disp(' ');
        
        % Link to ref
        if ~(replab.platformIsOctave  || (~usejava('desktop')))
            disp(['   <a href="matlab: ', helpFunctionName, '(''-f'',''', fullName, ''')">Reference page for ', fullName, '</a>'])
        end
        disp(' ');
    else
        replab.infra.dispH([' ', upper(k(1)), k(2:end), ' ', elementName, ':'], fullName, helpFunctionName, fullMode);
        disp(' ');
        
        % Class
        disp('   In class:');
        replab.infra.dispH(['     ', class.fullName], fullName, helpFunctionName, fullMode);
        disp(' ');
        
        % Also present in
        if length(elements) > 1
            disp('   Also present in:');
            for i = 2:length(elements)
                eli = elements{i};
                replab.infra.dispH(['     ', eli.classFullName], fullName, helpFunctionName, fullMode);
            end
            disp(' ');
        end
        
        % Declaration
        if isequal(k, 'method')
            el1 = elements{1};
            replab.infra.dispH(['   Declaration in ' el1.classFullName ':'], fullName, helpFunctionName, fullMode);
            replab.infra.dispH(['     ', el1.declaration], fullName, helpFunctionName, fullMode);
            disp(' ');
        end
        
        % type
        if isequal(k, 'property') && ~isempty(doc)
            type = doc.doc.firstLine;
            semicols = strfind(type, ':');
            if isempty(semicols)
                type = '';
            else
                type = type(1:semicols(1)-1);
            end
            if ~isempty(type)
                disp(['   Type: ', type]);
                disp(' ');
            end
        end
        
        % Attributes
        if isequal(k, 'property')
            if ~isempty(doc) && ~isempty(doc.attributes) && isfield(doc.attributes, 'SetAccess')
                disp(['   Attributes: SetAccess ', doc.attributes.SetAccess]);
                disp(' ');
            end
        elseif isequal(k, 'method')
            % TODO : display method attributes
%             el1 = elements{1};
%             if ~isempty(el1.attributes)
%                 el1.attributes;
%             end
        end
        
        % Link to help
        if ~(replab.platformIsOctave  || (~usejava('desktop')))
            disp(['   <a href="matlab: ', helpFunctionName, '(''', fullName, ''')">Help page for ', fullName, '</a>'])
        end
        disp(' ');
        
    end
    
    % Link to source code
    disp(sprintf('   <a href="matlab: opentoline(''%s'',%d,1)">See source</a>', class.fullFilename, elements{1}.lineNumber));
    disp(' ');
end

function help_function(codeBase, fun, helpFunctionName, fullMode)
    fullName = [strjoin(fun.packageNameParts, '.'), '.', fun.name];
    if ~fullMode
        replab.infra.dispH([' ', fun.declaration], fullName, helpFunctionName, fullMode);
        
        % Doc
        fun.doc.dispFilteredLines(fullName, helpFunctionName, fullMode);
        disp(' ');
        
        % Link to ref
        if ~(replab.platformIsOctave  || (~usejava('desktop')))
            disp(['   <a href="matlab: ', helpFunctionName, '(''-f'',''', fullName, ''')">Reference page for ', fullName, '</a>'])
        end
        disp(' ');
    else
        replab.infra.dispH([' Function ', fun.name, ':'], fullName, helpFunctionName, fullMode);
        disp(' ');
        
        % Package:
        if ~isempty(fun.packageNameParts)
            disp('   In package:');
            packageName = strjoin(fun.packageNameParts, '.');
            if replab.platformIsOctave  || (~usejava('desktop'))
                disp(['     ', packageName]);
            else
                disp(['     <a href="matlab: ', helpFunctionName, '(''-f'',''', packageName, ''')">', packageName, '</a>']);
            end
            disp(' ');
        end
        
        % Link to help
        if ~(replab.platformIsOctave  || (~usejava('desktop')))
            disp(['   <a href="matlab: ', helpFunctionName, '(''', fullName, ''')">Help page for ', fullName, '</a>'])
        end
        disp(' ');
    end
    
    % Link to source code
    disp(sprintf('   <a href="matlab: opentoline(''%s'',%d,1)">See source</a>', fun.fullFilename, 1));
    disp(' ');
end
