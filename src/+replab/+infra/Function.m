classdef Function < replab.infra.PackageElement
% Describes a MATLAB function
%
% Is also used to parse methods, thus the `parseAbstractBody` parsing node
    
    properties
        declaration
    end
    
    methods
        
        function self = Function(name, declaration, doc, packageNameParts)
            self.name = name;
            self.declaration = declaration;
            self.doc = doc;
            self.packageNameParts = packageNameParts;
            self.kind = 'function';
        end
        
    end
    
    methods (Static)
        
        function name = nameFromDeclaration(declaration)
        % Retrievs the function name from its declaration line
        %
        % Args:
        %   declaration (charstring): Trimmed function/method declaration line
            if sum(declaration == '=') >= 1
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
        
        function res = parseBlank(ct, pos)
            res = ct.expect(pos, ' ');
        end
        
        function res = parseComment(ct, pos)
            res = ct.expect(pos, '%');
        end
        
        function res = parseCode(ct, pos)
            res = ct.expect(pos, '!');
        end
        
        function pos = parseAbstractBody(ct, pos)
        % Parses the body of an abstract method, used for methods and not functions
            [pos line] = ct.expect(pos, '!');
            if isempty(pos) || ~isequal(strtrim(line), 'error(''Abstract'');')
                pos = [];
                return
            end
            pos = ct.expect(pos, '<');
        end
        
        function pos = parseControlStructure(ct, pos)
        % Parses a control structure such as 'if' or 'while'
        %
        % Assumes that the first line has already been consumed, and consumes the
        % rest of a control structure, including the final 'end'
            while 1
                res = replab.infra.Function.parseFunctionElement(ct, pos);
                if isempty(res)
                    break
                else
                    pos = res;
                end
            end
            pos = ct.expect(pos, '<');
        end
        
        function pos = parseFunctionElement(ct, pos)
        % Parses an element appearing in a function or method body
            tag = ct.peek(pos);
            switch tag
              case '$'
                error('End of file encountered: functions must end with an end statement');
              case {'c', 'p', 'm', 'e'}
                error(sprintf('classdef/properties/methods or enumeration are all invalid inside a function or method'));
              case {' ', '%', '!'}
                pos = ct.take(pos);
              case {'f', '>'}
                pos = ct.take(pos); % consume token
                pos = replab.infra.Function.parseControlStructure(ct, pos);
              otherwise
                pos = [];
            end
        end
        
        function [pos name declaration docLines isAbstract] = parse(ct, pos)
        % Parses a function or a method and returns the information fields
            name = [];
            declaration = [];
            docLines = {};
            isAbstract = false;
            [pos declaration] = ct.expect(pos, 'f');
            if isempty(pos)
                return
            end
            assert(~isempty(pos));
            name = replab.infra.Function.nameFromDeclaration(declaration);
            [pos docLines] = replab.infra.parseDocLines(ct, pos);
            res = replab.infra.Function.parseAbstractBody(ct, pos);
            if ~isempty(res)
                isAbstract = true;
                pos = res;
            else
                pos = replab.infra.Function.parseControlStructure(ct, pos);
            end
        end
        
        function f = fromParseState(ct, packageNameParts)
        % Parses a function and returns a `replab.infra.Function` instance
            pos = 1;
            [pos name declaration docLines isAbstract] = replab.infra.Function.parse(ct, pos);
            assert(~isempty(pos));
            assert(~isAbstract);
            doc = replab.infra.Doc.leftTrimmed(docLines);
            f = replab.infra.Function(name, declaration, doc, packageNameParts);
        end
        
    end
    
end
