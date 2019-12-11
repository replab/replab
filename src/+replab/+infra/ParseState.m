classdef ParseState < replab.Str
% Stores the contents and the parsing position (=state) of a source code file
%
% We process the source code files line by line, i.e. each line is a token
% for the parser; whitespace is trimmed from the lines we read.
%
% Each line has an associated tag.
% 
% Possible tags are:
%
% - 'BLANK': The line contains only whitespace
% - 'COMMENT': The line contains only a comment
% - 'CODE': The line contains code
% - 'CLASSDEF': The line contains a classdef statement
% - 'PROPERTIES': The line contains a methods statement
% - 'METHODS': The line contains a methods statement
% - 'ENUMERATION': The line contains an enumeration statement
% - 'FUNCTION': The line contains a function statement
% - 'IF': The line contains an if statement
% - 'TRY': The line contains a try statement
% - 'WHILE': The line contains a while statement
% - 'SWITCH': The line contains a switch statement
% - 'PARFOR': The line contains a parfor statement
% - 'SPMD': The line contains a spmd statement
% - 'FOR': The line contains a for statement
% - 'END': The line contains an end statement
% - 'EOF': The end of file token

    properties
        lines % row cell array of charstring: Source code lines
        tags % row cell array of charstring: Tag describing the line type
        pos % integer: Current line position
    end
    
    methods
        
        function self = ParseState(lines, tags, pos)
            self.tags = tags;
            self.lines = lines;
            self.pos = pos;
        end
        
        function [nextParseState tag line] = take(self)
        % Consumes a line from the input, and returns the updated parser state
        %
        % `tag` and `line` correspond to the consumed line.
            tag = self.tags{self.pos};
            if ~isequal(tag, 'EOF')
                line = self.lines{self.pos};
            else
                line = [];
            end
            nextParseState = replab.infra.ParseState(self.lines, self.tags, self.pos + 1);
        end
        
        function [nextParseState line] = expect(self, expectedTag)
        % Equivalent to `take`, conditioned on the consumed tag to be `expectedTag`
            [nextParseState tag line] = self.take;
            if ~isequal(tag, expectedTag)
                nextParseState = [];
                line = [];
            end
        end
        
        function tag = peek(self)
        % Peeks at the tag of the next line to be consumed, without consuming it
            tag = self.tags{self.pos};
        end
        
    end
    
    methods (Static)
       
        function ps = fromSourceLines(lines)
        % Constructs a ParseState instance from source code lines
        %
        % Args:
        %   lines (row cell array of charstring): Trimmed source code lines
        %
        % Returns:
        %   :class:`replab.infra.ParseState`: A fresh ParseState instance
            n = length(lines);
            tags = cell(1, n+1);
            tags{n+1} = 'EOF';
            for i = 1:n
                line = lines{i};
                if isempty(line)
                    tags{i} = 'BLANK';
                elseif line(1) == '%'
                    tags{i} = 'COMMENT';
                else
                    word = regexp(line, '^\w+', 'match');
                    if isempty(word)
                        tags{i} = 'CODE';
                    else
                        switch word{1}
                          case {'classdef', 'properties', 'methods', 'enumeration', 'function', 'if', 'try', 'while', ...
                                'switch', 'parfor', 'spmd', 'for', 'end'}
                            tags{i} = upper(word{1});
                          otherwise
                            tags{i} = 'CODE';
                        end
                    end
                end
            end
            ps = replab.infra.ParseState(lines, tags, 1);
        end
        
        function ps = fromSource(contents)
        % Constructs a ParseState instance from the source code
        %
        % Args:
        %   contents (charstring): File contents as a single char string
        %
        % Returns:
        %   :class:`replab.infra.ParseState`: A fresh ParseState instance
            lines = cellfun(@strtrim, regexp(contents, '[\n\r]+', 'split'), 'uniform', 0);
            ps = replab.infra.ParseState.fromSourceLines(lines);
        end
        
        function ps = fromFile(filename)
        % Constructs a ParseState instance from the contents of the given filename
        %
        % Args:
        %   filename (charstring): Name of the source file to read and lex
        %
        % Returns:
        %   :class:`replab.infra.ParseState`: A fresh ParseState instance
            contents = fileread(filename);
            ps = replab.infra.ParseState.fromSource(contents);
        end
        
    end
    
end
