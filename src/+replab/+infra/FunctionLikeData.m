classdef FunctionLikeData < replab.Str
% Describes the data recovered when parsing a MATLAB function or method
%
% The code below also contains the parser static methods to parse that data, and
% that parsing code is also used for class methods; thus the presence of a `parseAbstractBody` 
% parsing node, and a isAbstract property.
    
    properties
        name % charstring: Function or method name
        declaration % charstring: Function or method declaration line
        declarationLineNumber % integer: Line number of the declaration
        docLines % row cell array of charstring: Documentation comment lines stripped of leading whitespace and leading ``%``
        docLineNumbers % row integer vector: Line numbers of the documentation comment
        attributes % struct: Attributes from the ``methods`` block (or ``[]`` if this describes a function)
    end
    
    methods
        
        function self = FunctionLikeData(name, declaration, declarationLineNumber, ...
                                         docLines, docLineNumbers, attributes)
            self.name = name;
            self.declaration = declaration;
            self.declarationLineNumber = declarationLineNumber;
            self.docLines = docLines;
            self.docLineNumbers = docLineNumbers;
            self.attributes = attributes;
        end
        
    end
    
    methods (Static)
        
        function name = nameFromDeclaration(ct, pos, declaration)
        % Retrievs the function name from its declaration line
        %
        % Args:
        %   declaration (charstring): Trimmed function/method declaration line
            if sum(declaration == '=') >= 1
                parts = strsplit(declaration, '=');
                if length(parts) ~= 2
                    replab.infra.parseError(ct, pos, 'Invalid number of =');
                end
                tokens = regexp(strtrim(parts{2}), '^(\w+)', 'tokens', 'once');
                if length(tokens) ~= 1
                    replab.infra.parseError(ct, pos, 'Cannot find method/function name in declaration');
                end
                name = tokens{1};
            else
                tokens = regexp(declaration, '^function\s+(\w+)', 'tokens', 'once');
                if length(tokens) ~= 1
                    replab.infra.parseError(ct, pos, 'Cannot find method/function name in declaration');
                end                
                name = tokens{1};
            end
            if ~isempty(regexp(name, '__'))
                replab.infra.parseError(ct, pos, 'Name of method/function cannot contain double underscore');
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
                res = replab.infra.FunctionLikeData.parseFunctionElement(ct, pos);
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
                replab.infra.parseError(ct, pos, 'Functions must end with an end statement');
              case {'c', 'p', 'm', 'e'}
                replab.infra.parseError(ct, pos, 'Invalid token inside a function or method');
              case {' ', '%', '!'}
                pos = ct.take(pos);
              case {'f', '>'}
                pos = ct.take(pos); % consume token
                pos = replab.infra.FunctionLikeData.parseControlStructure(ct, pos);
              otherwise
                pos = [];
            end
        end
        
        function [pos fld] = parse(ct, pos, attributes)
        % Parses a function or a method and returns the information
            fld = [];
            startPos = pos;
            [pos declaration] = ct.expect(pos, 'f');
            if isempty(pos)
                return
            end
            name = replab.infra.FunctionLikeData.nameFromDeclaration(ct, startPos, declaration);
            [pos docLines docLineNumbers] = replab.infra.parseDocLines(ct, pos);
            res = replab.infra.FunctionLikeData.parseAbstractBody(ct, pos);
            if ~isempty(res)
                % The function/method has an abstract body
                if isempty(attributes)
                    replab.infra.parseError(ct, pos, 'Only methods can have an abstract body');
                else
                    attributes.Abstract = true;
                end
                pos = res;
            else
                % The function/method does not have an abstract body, so parse a real body instead
                pos = replab.infra.FunctionLikeData.parseControlStructure(ct, pos);
            end
            fld = replab.infra.FunctionLikeData(name, declaration, startPos, docLines, docLineNumbers, attributes);
        end
        
    end
    
end
