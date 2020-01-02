classdef ParseState < replab.Str
% Stores the contents and the parsing position (=state) of a doctest block
%
% We process the doctest block line by line, i.e. each line is a token for the parser;
% whitespace is trimmed from the lines we read, and empty lines are removed.
%
% For source code lines ('START' or 'CONT'), we also split the line into a code part and a comment part.
%
% Each line has an associated tag.
%
% Possible tags are:
%
% - 'START': The line starts with >>>
% - 'CONT': The line starts with ...
% - 'OUT': The line starts with another prefix
% - 'EOF': End of the block
%
% See https://jayconrod.com/posts/65/how-to-build-a-parser-by-hand
% for another example of a hand-written parser

    properties
        tags % row cell array of charstring: Tag describing the line type
        lines % row cell array of charstring: Content lines
        comments % row cell array of charstring: Comments
        lineNumbers % row vector of integer: Original position of remaining lines
        pos % integer: Current line position
    end

    methods

        function self = ParseState(tags, lines, comments, lineNumbers, pos)
            self.tags = tags;
            self.lines = lines;
            self.lineNumbers = lineNumbers;
            self.comments = comments;
            self.pos = pos;
        end

        function [nextParseState tag line comment lineNumber] = take(self)
        % Consumes a line from the input, and returns the updated parser state
        %
        % `tag` and `line` correspond to the consumed line.
            tag = self.tags{self.pos};
            line = self.lines{self.pos};
            comment = self.comments{self.pos};
            lineNumber = self.lineNumbers(self.pos);
            nextParseState = replab.infra.doctests.ParseState(self.tags, self.lines, self.comments, self.lineNumbers, self.pos + 1);
        end

        function [nextParseState line comment lineNumber] = expect(self, expectedTag)
        % Equivalent to `take`, conditioned on the consumed tag to be ``expectedTag``
            [nextParseState tag line comment lineNumber] = self.take;
            if ~isequal(tag, expectedTag)
                nextParseState = [];
            end
        end

    end

    methods (Static)

        function ps = fromDocTestBlock(inputLines)
        % Constructs a ParseState instance from a doctest block lines
        %
        % Args:
        %   lines (row cell array of charstring): Doctest block lines
        %
        % Returns:
        %   `.ParseState`: A fresh ParseState instance
            lines = {};
            comments = {};
            tags = {};
            lineNumbers = [];
            for i = 1:length(inputLines)
                l = strtrim(inputLines{i});
                if isempty(l) || l(1) == '%'
                    % omit that line
                elseif length(l) >= 3 && isequal(l(1:3), '>>>')
                    tags{1,end+1} = 'START';
                    tmp = l(4:end);
                    [code comment] = replab.infra.doctests.splitComment(tmp);
                    lines{1,end+1} = strtrim(code);
                    comments{1,end+1} = strtrim(comment);
                    lineNumbers(1,end+1) = i;
                elseif length(l) >= 3 && isequal(l(1:3), '...')
                    tags{1,end+1} = 'CONT';
                    tmp = l(4:end);
                    [code comment] = replab.infra.doctests.splitComment(tmp);
                    lines{1,end+1} = strtrim(code);
                    comments{1,end+1} = strtrim(comment);
                    lineNumbers(1,end+1) = i;
                else
                    tags{1,end+1} = 'OUT';
                    lines{1,end+1} = l;
                    comments{1,end+1} = '';
                    lineNumbers(1,end+1) = i;
                end
            end
            tags{1,end+1} = 'EOF';
            lines{1,end+1} = '';
            comments{1,end+1} = '';
            lineNumbers(1,end+1) = 0;
            ps = replab.infra.doctests.ParseState(tags, lines, comments, lineNumbers, 1);
        end

    end

end
