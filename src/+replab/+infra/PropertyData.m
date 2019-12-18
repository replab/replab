classdef PropertyData < replab.Str
% Describes the data recovered when parsing a MATLAB class property
    
    properties
        name % charstring: Property name
        declarationLineNumber % integer: Line number of the property declaration
        docLines % row cell array of charstring: Documentation comment lines stripped of leading whitespace and leading ``%``
        docLineNumbers % row integer vector: Line numbers of the documentation comment
                       %
                       %                     Overlaps with `declarationLineNumber`
        attributes % struct: Attributes from the ``properties`` block
    end
    
    methods

        function self = PropertyData(name, declarationLineNumber, docLines, docLineNumbers, attributes)
            self.name = name;
            self.declarationLineNumber = declarationLineNumber;
            self.docLines = docLines;
            self.docLineNumbers = docLineNumbers;
            self.attributes = attributes;
        end
        
    end
    
    methods (Static)
       
        function [pos pd] = parse(ct, pos, attributes)
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
            startPos = pos;
            [pos line] = ct.expect(pos, '!');
            if isempty(pos)
                pd = [];
                return
            end
            
            % splits the property line around a possible % indicating a comment
            parts = strsplit(line, '%');
            def = parts{1};
            firstDocLine = strjoin(parts(2:end), '%');
            
            % splits the property code to get the property name (i.e. anything before the first = or ;)
            parts = regexp(def, '[=;]', 'split');
            name = strtrim(parts{1});
            if isempty(name)
                replab.infra.parseError(ct, startPos, 'Cannot find property name');
            end
            underpos = regexp(name, '_[^_]');
            if ~isempty(underpos)
                replab.infra.parseError(ct, pos, 'Name of property cannot contain underscore except in terminal position');
            end
            if isempty(firstDocLine)
                docLines = {};
                docLineNumbers = [];
            else
                docLines = {firstDocLine};
                docLineNumbers = startPos;
                l = 2;
                while 1
                    [res line] = ct.expect(pos, '%');
                    if isempty(res)
                        break
                    else
                        content = line(2:end);
                        if l == 2 && ~isempty(strtrim(content))
                            replab.infra.parseWarning(ct, pos, 'Second documentation comment line should be empty');
                        end
                        docLines{1,l} = content;
                        docLineNumbers(1,l) = pos;
                        pos = res;
                    end
                    l = l + 1;
                end
            end
            pd = replab.infra.PropertyData(name, startPos, docLines, docLineNumbers, attributes);
        end
        
    end

end
