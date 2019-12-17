classdef ClassData < replab.Str

    properties
        name % charstring: Class name
        parentNames % row cell array of charstring: Full qualified names of parents
        docLines % row cell array of charstring: Documentation comment lines
                 %                               stripped of leading whitespace and leading ``%``
        docLineNumbers % row integer vector: Line numbers of the documentation comment
        members % row
    end
    
    methods
        
        function self = Class(name, parentNames, docLines, docLineNumbers, members)
            self.name = name;
            self.parentsNames = parentsNames;
            self.doc = doc;
            self.packageNameParts = packageNameParts;
            members = struct;
            for i = 1:length(memberList)
                member = memberList{i};
                members.(member.name) = member;
            end
            members = orderfields(members);
            self.members = members;
            self.kind = 'class';
            self.fullFilename = fullFilename;
        end
        
        function str = fullName(self)
            str = strjoin(horzcat(self.packageNameParts, {self.name}), '.');
        end

        function n = nParents(self)
            n = length(self.parentsNames);
        end
        
        function nm = parentName(self, i)
            nm = self.parentsNames{i};
        end
        
        function nm = parentNameParts(self, i)
            nm = strsplit(self.parentsNames{i}, '.');
        end
        
        function mn = memberNames(self)
            mn = fieldnames(self.members)';
        end
        
        function p = parent(self, codeBase, i)
            if isequal(self.parentName(i), 'handle')
                p = [];
            else
                nameParts = strsplit(self.parentName(i), '.');
                [package packageNameParts restNameParts] = codeBase.lookupGreedy(nameParts);
                assert(length(restNameParts) == 1);
                p = package.member(restNameParts{1});
            end
        end
        
        function c = childrenNames(self, codeBase)
        % Returns the names of the children of the current class
        % 
        % Args:
        %   codeBase (`.CodeBase`): Code base containing this class
        %
        % Returns:
        %   row cell vector of charstring: Names of the subclasses of the current class
            fne = self.fieldNameEncoding;
            if isfield(codeBase.subclasses, fne)
                c = codeBase.subclasses.(fne);
            else
                c = {};
            end
        end
        
        function b = hasMember(self, name)
            b = isfield(self.members, name);
        end
        
        function b = hasInheritedMember(self, codeBase, name)
            b = true;
            if self.hasMember(name)
                return
            end
            for i = 1:self.nParents
                p = self.parent(codeBase, i);
                if ~isempty(p)
                    if p.hasInheritedMember(codeBase, name)
                        return
                    end
                end
            end
            b = false;
        end
        
        function m = member(self, name)
            m = self.members.(name);
        end
        
        function str = fieldNameEncoding(self)
            str = strjoin(horzcat(self.packageNameParts, {self.name}), '_');
        end
        
        function elements = findInheritedMember(self, codeBase, name)
            elements = {};
            visitedClassIds = {};
            stack = {self};
            % breath-first search through the inheritance tree
            while ~isempty(stack)
                % pop item off the stack
                h = stack{1};
                stack = stack(2:end);
                id = h.fieldNameEncoding;
                if h.hasMember(name) && ~ismember(id, visitedClassIds)
                    elements{end+1} = h.member(name);
                    visitedClassIds{end+1} = id;
                end
                for i = 1:h.nParents
                    p = h.parent(codeBase, i);
                    if ~isempty(p)
                        stack{end+1} = p;
                    end
                end
            end
        end
        
        function names = inheritedMemberNames(self, codeBase)
            names = {};
            stack = {self};
            % breath-first search through the inheritance tree
            while ~isempty(stack)
                % pop item off the stack
                h = stack{1};
                stack = stack(2:end);
                % add names of popped item members
                names = union(names, h.memberNames);
                % add parents to the stack
                for i = 1:h.nParents
                    p = h.parent(codeBase, i);
                    if ~isempty(p)
                        stack{end+1} = p;
                    end
                end
            end
        end
        
        function names = hiddenFields(self)
            names = hiddenFields@replab.Str(self);
            names{1, end+1} = 'members';
        end
        
        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Str(self);
            fn = fieldnames(self.members);
            for i = 1:length(fn)
                names{1, end+1} = sprintf('member(''%s'')', fn{i});
                values{1, end+1} = self.member(fn{i});
            end
        end

    end
    
    methods (Static)
        
        function [pos members] = parseMethod(ct, pos, packageNameParts, className, attributes)
        % Parses a method declaration
            [pos fld] = replab.infra.FunctionLikeData.parse(ct, pos, attributes);
            if isempty(pos)
                members = {};
                return
            end
            doc = replab.infra.Doc.leftTrimmed(fld.docLines);
            members = {replab.infra.Method(fld.name, attributes, fld.declaration, doc, ...
                                           packageNameParts, className, fld.declarationLineNumber)};
        end

        function [res members] = parseMethodsElement(ct, pos, packageNameParts, className, attributes)
        % Parses an element that can appear in a methods block
            [res members] = replab.infra.Class.parseBlankAsEmptyCell(ct, pos);
            if ~isempty(res)
                return
            end
            [res members] = replab.infra.Class.parseCommentAsEmptyCell(ct, pos);
            if ~isempty(res)
                return
            end
            [res members] = replab.infra.Class.parseMethod(ct, pos, packageNameParts, className, attributes);
            if ~isempty(res)
                return
            end
            res = [];
            members = {};
        end

        function [pos members] = parseMethods(ct, pos, packageNameParts, className)
        % Parses a methods block
            members = {};
            [pos line] = ct.expect(pos, 'm');
            if isempty(pos)
                return
            end
            % Remove an eventual comment
            parts = strsplit(line, '%');
            line = strtrim(parts{1});
            if isequal(line, 'methods')
                % No attributes
                attributes = struct;
            else
                % Has attributes
                tokens = regexp(line, '^methods\s*\((.*)\)$', 'tokens', 'once');
                assert(length(tokens) == 1);
                attributes = replab.infra.parseAttributes(tokens{1});
            end
            while 1
                [res newMembers] = replab.infra.Class.parseMethodsElement(ct, pos, packageNameParts, className, attributes);
                if isempty(res)
                    break
                else
                    pos = res;
                    members = horzcat(members, newMembers);
                end
            end
            [pos line] = ct.expect(pos, '<');
            assert(~isempty(pos));
        end

        function [pos members] = parseProperty(ct, pos, attributes, packageNameParts, className, filename)
            [pos pd] = replab.infra.PropertyData.parse(ct, pos, attributes);
            if isempty(pos)
                members = {};
                return
            end
            doc = replab.infra.Doc.leftTrimmed(pd.docLines);
            property = replab.infra.Property(pd.name, pd.attributes, doc, packageNameParts, className, pd.declarationLineNumber);
            members = {property};
        end
        
        function [res members] = parsePropertiesElement(ct, pos, attributes, packageNameParts, className)
        % Parses an element that can appear in a properties block
            [res members] = replab.infra.Class.parseBlankAsEmptyCell(ct, pos);
            if ~isempty(res)
                return
            end
            [res members] = replab.infra.Class.parseCommentAsEmptyCell(ct, pos);
            if ~isempty(res)
                return
            end
            [res members] = replab.infra.Class.parseProperty(ct, pos, attributes, packageNameParts, className);
            if ~isempty(res)
                return
            end
            res = [];
            members = {};
        end

        function [pos members] = parseProperties(ct, pos, packageNameParts, className)
        % Parses a properties block
            members = {};
            [pos line] = ct.expect(pos, 'p');
            if isempty(pos)
                return
            end
            % Remove an eventual comment
            parts = strsplit(line, '%');
            line = strtrim(parts{1});
            if isequal(line, 'properties')
                % No attributes
                attributes = struct;
            else
                % Has attributes
                tokens = regexp(line, '^properties\s*\((.*)\)$', 'tokens', 'once');
                assert(length(tokens) == 1);
                attributes = replab.infra.parseAttributes(tokens{1});
            end
            while 1
                [res newMembers] = replab.infra.Class.parsePropertiesElement(ct, pos, attributes, packageNameParts, className);
                if isempty(res)
                    break
                else
                    pos = res;
                    members = horzcat(members, newMembers);
                end
            end
            [pos line] = ct.expect(pos, '<');
            assert(~isempty(pos));
        end
        
        function [res newMembers] = parseBlankAsEmptyCell(ct, pos)
        % Parses a blank line returning an empty list of members
            res = ct.expect(pos, ' ');
            newMembers = {};
        end

        function [res newMembers] = parseCommentAsEmptyCell(ct, pos)
        % Parses a comment line returning an empty list of members
            res = ct.expect(pos, '%');
            newMembers = {};
        end
        
        function [res members] = parseClassElement(ct, pos, packageNameParts, className)
        % Parses an element that can appear in a class definition
            [res members] = replab.infra.Class.parseBlankAsEmptyCell(ct, pos);
            if ~isempty(res)
                return
            end
            [res members] = replab.infra.Class.parseCommentAsEmptyCell(ct, pos);
            if ~isempty(res)
                return
            end
            [res members] = replab.infra.Class.parseProperties(ct, pos, packageNameParts, className);
            if ~isempty(res)
                return
            end
            [res members] = replab.infra.Class.parseMethods(ct, pos, packageNameParts, className);
            if ~isempty(res)
                return
            end
            res = [];
            members = {};
        end
        
        function [className parentsNames docLines memberList] = parse(ct, packageNameParts)
        % Parses a full class definition
            pos = 1;
            [pos line] = ct.expect(pos, 'c');
            assert(~isempty(pos));
            tokens = cellfun(@strtrim, regexp(line, '[&<]', 'split'), 'uniform', 0);
            matches = regexp(tokens{1}, '^classdef\s+(\w+)$', 'tokens', 'once');
            className = matches{1};
            parentsNames = tokens(2:end);
            [pos docLines] = replab.infra.parseDocLines(ct, pos);
            assert(~isempty(pos));
            memberList = {};
            while 1
                [res newMembers] = replab.infra.Class.parseClassElement(ct, pos, packageNameParts, className);
                if isempty(res)
                    break
                else
                    pos = res;
                    memberList = horzcat(memberList, newMembers);
                end
            end
            pos = ct.expect(pos, '<');
            assert(~isempty(pos));
        end
        
        function c = fromParseState(ct, packageNameParts, filename)
        % Parses a full definition and constructs a `+replab.+infra.Class` instance
            [className parentsNames docLines memberList] = replab.infra.Class.parse(ct, packageNameParts);
            doc = replab.infra.Doc.leftTrimmed(docLines);
            c = replab.infra.Class(className, parentsNames, doc, memberList, packageNameParts, filename);
        end
        
    end
    
end
