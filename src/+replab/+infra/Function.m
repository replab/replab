classdef Function < replab.infra.PackageElement
% Describes a MATLAB function
%
% Is also used to parse methods, thus the `parseAbstractBody` parsing node
    
    properties
        declaration
    end
    
    methods
        
        function self = Function(name, declaration, docLines)
            self.name = name;
            self.declaration = declaration;
            self.docLines = docLines;
        end
        
        function str = headerStr(self)
            str = sprintf('%s (function)', self.name);
        end
        
    end
    
    methods (Static)
        
        function name = nameFromDeclaration(declaration)
        % Retrievs the function name from its declaration line
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
        
        function ps = parseAbstractBody(ps)
        % Parses the body of an abstract method, used for methods and not functions
            [ps line] = ps.expect('CODE');
            if isempty(ps) || ~isequal(strtrim(line), 'error(''Abstract'');')
                ps = [];
                return
            end
            ps = ps.expect('END');
        end
        
        function ps = parseControlStructure(ps)
        % Parses a control structure such as 'if' or 'while'
        %
        % Assumes that the first line has already been consumed, and consumes the
        % rest of a control structure, including the final 'end'
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
        % Parses an element appearing in a function or method body
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
        
        function [ps name declaration docLines isAbstract] = parse(ps)
        % Parses a function or a method and returns the information fields
            name = [];
            declaration = [];
            docLines = {};
            isAbstract = false;
            [ps declaration] = ps.expect('FUNCTION');
            if isempty(ps)
                return
            end
            assert(~isempty(ps));
            name = replab.infra.Function.nameFromDeclaration(declaration);
            [ps docLines] = replab.infra.parseDocLines(ps);
            res = replab.infra.Function.parseAbstractBody(ps);
            if ~isempty(res)
                isAbstract = true;
                ps = res;
            else
                ps = replab.infra.Function.parseControlStructure(ps);
            end
        end
        
        function f = fromParseState(ps)
        % Parses a function and returns a `replab.infra.Function` instance
            [ps name declaration docLines isAbstract] = replab.infra.Function.parse(ps);
            assert(~isempty(ps));
            assert(~isAbstract);
            f = replab.infra.Function(name, declaration, docLines);
        end
        
    end
    
end
