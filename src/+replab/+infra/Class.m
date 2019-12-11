classdef Class < replab.Str

    properties
        className
        parentsNames
        docLines
        members
        %myMethods
        %myProperties
    end
    
    methods
        
        function self = Class(className, parentsNames, docLines, members)
            self.className = className;
            self.parentsNames = parentsNames;
            self.docLines = docLines;
            self.members = members;
        end
        
        function n = nMembers(self)
            n = length(self.members);
        end
        
        function m = member(self, i)
            m = self.members{i};
        end
        
        function names = hiddenFields(self)
            names = hiddenFields@replab.Str(self);
            names{1, end+1} = 'members';
        end
        
        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Str(self);
            for i = 1:self.nMembers
                names{1, end+1} = sprintf('member(%d)', i);
                values{1, end+1} = self.member(i);
            end
        end
        
        function member = lookupMemberName(self, name)
            
        end
        
    end
    
    methods (Static)

        function [ps members] = parseMethod(ps, attributes)
            [ps name declaration docLines isAbstract] = replab.infra.Function.parse(ps);
            if isempty(ps)
                members = {};
                return
            end
            members = {replab.infra.Method(name, attributes, declaration, docLines, isAbstract)};
        end

        function [res members] = parseMethodsElement(ps, attributes)
            [res members] = replab.infra.Class.parseBlankAsEmptyCell(ps);
            if ~isempty(res)
                return
            end
            [res members] = replab.infra.Class.parseCommentAsEmptyCell(ps);
            if ~isempty(res)
                return
            end
            [res members] = replab.infra.Class.parseMethod(ps, attributes);
            if ~isempty(res)
                return
            end
            res = [];
            members = {};
        end

        function [ps members] = parseMethods(ps)
            members = {};
            [ps line] = ps.expect('METHODS');
            if isempty(ps)
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
                [res newMembers] = replab.infra.Class.parseMethodsElement(ps, attributes);
                if isempty(res)
                    break
                else
                    ps = res;
                    members = horzcat(members, newMembers);
                end
            end
            [ps line] = ps.expect('END');
            assert(~isempty(ps));
        end

        function [ps members] = parseProperty(ps, attributes)
        % Parses a property definition
        %
        % The formats can be
        %
        % name
        % name;
        % name = value
        % name = value;
        %
        % followed eventually by a comment as in
        %
        % name = value; % documentation comment
        %
        % and possibly multiline comments as in
        %
        % name = value; % documentation comment
        %               % comment continued
            [ps line] = ps.expect('CODE');
            if isempty(ps)
                members = {};
                return
            end
            
            % splits the property line around a possible % indicating a comment
            parts = strsplit(line, '%');
            def = parts{1};
            firstDocLine = strjoin(parts(2:end), '%');
            
            % splits the property code to get the property name (i.e. anything before the first = or ;)
            parts = regexp(def, '[=;]', 'split');
            name = strtrim(parts{1});
            assert(~isempty(name));
            if isempty(firstDocLine)
                docLines = {};
            else
                nextDocLines = {};
                while 1
                    [res line] = ps.expect('COMMENT');
                    if isempty(res)
                        break
                    else
                        nextDocLines{1, end+1} = line(2:end);
                        ps = res;
                    end
                end
                docLines = horzcat({firstDocLine}, nextDocLines);
            end
            property = replab.infra.Property(name, attributes, docLines);
            members = {property};
        end
        
        function [res members] = parsePropertiesElement(ps, attributes)
            [res members] = replab.infra.Class.parseBlankAsEmptyCell(ps);
            if ~isempty(res)
                return
            end
            [res members] = replab.infra.Class.parseCommentAsEmptyCell(ps);
            if ~isempty(res)
                return
            end
            [res members] = replab.infra.Class.parseProperty(ps, attributes);
            if ~isempty(res)
                return
            end
            res = [];
            members = {};
        end

        function [ps members] = parseProperties(ps)
            members = {};
            [ps line] = ps.expect('PROPERTIES');
            if isempty(ps)
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
                [res newMembers] = replab.infra.Class.parsePropertiesElement(ps, attributes);
                if isempty(res)
                    break
                else
                    ps = res;
                    members = horzcat(members, newMembers);
                end
            end
            [ps line] = ps.expect('END');
            assert(~isempty(ps));
        end
        
        function [res newMembers] = parseBlankAsEmptyCell(ps)
            res = ps.expect('BLANK');
            newMembers = {};
        end

        function [res newMembers] = parseCommentAsEmptyCell(ps)
            res = ps.expect('COMMENT');
            newMembers = {};
        end
        
        function [res members] = parseClassElement(ps)
            [res members] = replab.infra.Class.parseBlankAsEmptyCell(ps);
            if ~isempty(res)
                return
            end
            [res members] = replab.infra.Class.parseCommentAsEmptyCell(ps);
            if ~isempty(res)
                return
            end
            [res members] = replab.infra.Class.parseProperties(ps);
            if ~isempty(res)
                return
            end
            [res members] = replab.infra.Class.parseMethods(ps);
            if ~isempty(res)
                return
            end
            res = [];
            members = {};
        end
        
        function [className parentsNames docLines members] = parse(ps)
            [ps line] = ps.expect('CLASSDEF');
            assert(~isempty(ps));
            tokens = cellfun(@strtrim, regexp(line, '[&<]', 'split'), 'uniform', 0);
            matches = regexp(tokens{1}, '^classdef\s+(\w+)$', 'tokens', 'once');
            className = matches{1};
            parentsNames = tokens(2:end);
            [ps docLines] = replab.infra.parseDocLines(ps);
            assert(~isempty(ps));
            members = {};
            while 1
                [res newMembers] = replab.infra.Class.parseClassElement(ps);
                if isempty(res)
                    break
                else
                    ps = res;
                    members = horzcat(members, newMembers);
                end
            end
            [ps line] = ps.expect('END');
            assert(~isempty(ps));
        end
        
        function c = fromParseState(ps)
            [className parentsNames docLines members] = replab.infra.Class.parse(ps);
            c = replab.infra.Class(className, parentsNames, docLines, members);
        end
        
    end
    
end
