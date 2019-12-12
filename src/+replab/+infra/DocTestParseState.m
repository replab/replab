classdef DocTestParseState < replab.Str
% Stores the contents and the parsing position (=state) of a doctest block 
%
% We process the doctest block line by line, i.e. each line is a token
% for the parser; whitespace is trimmed from the lines we read.
%
% Each line has an associated tag.
% 
% Possible tags are:
%
% - 'BLANK': The line contains only whitespace
% - 'START': The line starts with >>>
% - 'CONT': The line starts with ...
% - 'OUT': The line starts with another prefix
% - 'EOF': End of the block
%
% See https://jayconrod.com/posts/65/how-to-build-a-parser-by-hand
% for another example of a hand-written parser
    
    properties
        lines % row cell array of charstring: Source code lines
        tags % row cell array of charstring: Tag describing the line type
        pos % integer: Current line position
    end
    
    methods
        
        function self = DocTestParseState(lines, tags, pos)
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
       
        function ps = fromDocTestBlock(lines)
        % Constructs a ParseState instance from a doctest block lines
        %
        % Args:
        %   lines (row cell array of charstring): Trimmed source code lines
        %
        % Returns:
        %   :class:`replab.infra.ParseState`: A fresh ParseState instance
            lines = cellfun(@strtrim, lines, 'uniform', 0);
            n = length(lines);
            tags = cell(1, n+1);
            tags{n+1} = 'EOF';
            for i = 1:n
                line = lines{i};
                if isempty(line)
                    tags{i} = 'BLANK';
                elseif length(line) >= 3 && isequal(line(1:3), '>>>')
                    tags{i} = 'START';
                elseif length(line) >= 3 && isequal(line(1:3), '...')
                    tags{i} = 'CONT';
                else
                    tags{i} = 'OUT';
                end
            end
            ps = replab.infra.DocTestParseState(lines, tags, 1);
        end
        
    end
    
end
