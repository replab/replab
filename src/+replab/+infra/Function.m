classdef Function < replab.Str
    
    properties
        name
        declaration
        docLines
    end
   
    methods
        
        function self = Function(name, declaration, docLines)
            self.name = name;
            self.declaration = declaration;
            self.docLines = docLines;
        end
        
    end
    
    methods (Static)
        
        function name = nameFromDeclaration(declaration)
        %
        % Args:
        %   declaration (charstring): Trimmed function/method declaration line
            if contains(declaration, '=')
                parts = strsplit(declaration, '=');
                assert(length(parts) == 2);
                tokens = regexp(strtrim(parts{2}), '^(\w+)', 'tokens', 'once');
                assert(length(tokens) == 1);
                name = tokens{1};
            else
                tokens = regexp(declaration, '^function\s+(\w+)', 'tokens', 'once');
                name = tokens{1};
            end
        end
        
        function res = parseBlank(ps)
            res = ps.expect('BLANK');
        end
        
        function res = parseComment(ps)
            res = ps.expect('COMMENT');
        end
        
        function res = parseCode(ps)
            res = ps.expect('CODE');
        end
        
        function ps = parseControlStructure(ps)
        % Parses the rest of a control structure after the starting token has been consumed
        %
        % Consumes also the terminal 'end'
            while 1
                res = replab.infra.Function.parseFunctionElement(ps);
                if isempty(res)
                    break
                else
                    ps = res;
                end
            end
            ps = ps.expect('END');
        end
        
        function ps = parseFunctionElement(ps)
            tag = ps.peek;
            switch tag
              case 'EOF'
                error('End of file encountered: functions must end with an end statement');
              case {'CLASSDEF', 'PROPERTIES', 'METHODS', 'ENUMERATION'}
                error(sprintf('Token %s invalid inside a function or method', tag));
              case {'BLANK', 'COMMENT', 'CODE'}
                ps = ps.take;
              case {'FUNCTION', 'IF', 'TRY', 'WHILE', 'SWITCH', 'PARFOR', 'SPMD', 'FOR'}
                ps = ps.take; % consume token
                ps = replab.infra.Function.parseControlStructure(ps);
              otherwise
                ps = [];
            end
        end
        
        function [ps name declaration docLines] = parse(ps)
            name = [];
            declaration = [];
            docLines = {};
            [ps declaration] = ps.expect('FUNCTION');
            if isempty(ps)
                return
            end
            assert(~isempty(ps));
            name = replab.infra.Function.nameFromDeclaration(declaration);
            [ps docLines] = replab.infra.parseDocLines(ps);
            ps = replab.infra.Function.parseControlStructure(ps);
        end
        
        function f = fromParseState(ps)
            [ps name declaration docLines] = replab.infra.Function.parse(ps);
            assert(~isempty(ps));
            f = replab.infra.Function(name, declaration, docLines);
        end
        
    end
        
end
