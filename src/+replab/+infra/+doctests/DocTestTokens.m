classdef DocTestTokens < replab.Str
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
% - 's': The line starts with >>> and does not end with ...
% - 'S': The line starts with >>> and ends with ...
% - 'o': The line starts with another prefix and does not end with ...
% - 'O': The line starts with another prefix and ends with ...
% - '$': End of the block
%
% See https://jayconrod.com/posts/65/how-to-build-a-parser-by-hand for another example of a hand-written parser

    properties (SetAccess = protected)
        tags % (char(1,\*)): Tag describing the line type
        lines % (cell(1,\*) of charstring): Line content
        comments % (cell(1,\*) of charstring): Line comment
        lineNumbers % (integer(1,\*)): Original line number
    end

    methods

        function self = DocTestTokens(tags, lines, comments, lineNumbers)
            self.tags = tags;
            self.lines = lines;
            self.comments = comments;
            self.lineNumbers = lineNumbers;
        end

        function t = peek(self, pos)
        % Peeks at the tag of the next line to be consumed, without consuming it
        %
        % Args:
        %   pos (integer): Line position
        %
        % Returns:
        %   char: Tag
            t = self.tags(pos);
        end

        function [newPos, tag, line, comment, lineNumber] = take(self, pos)
        % Consumes a line from the input, and returns the updated parser state
        %
        % Returns
        % -------
        %   newPos: integer
        %     New line position
        %   tag: char
        %     Consumed line tag
        %   line: charstring
        %     Line contents, with comment stripped
        %   comment: charstring
        %     Comment; if no comment present, empty string
        %   lineNumber: integer
        %     Original line number
            newPos = pos + 1;
            tag = self.tags(pos);
            line = self.lines{pos};
            comment = self.comments{pos};
            lineNumber = self.lineNumbers(pos);
        end

    end

    methods (Static)

        function dtt = lex(inputLines)
        % Constructs a DocTestTokens instance from a doctest block lines
        %
        % Args:
        %   lines (row cell array of charstring): Doctest block lines
        %
        % Returns:
        %   `.ParseState`: A fresh ParseState instance
            lines = {};
            comments = {};
            tags = '';
            lineNumbers = [];
            ind = 1;
            for i = 1:length(inputLines)
                l = strtrim(inputLines{i});
                if isempty(l) || l(1) == '%'
                    % omit that line
                elseif replab.compat.startsWith(l, '>>>')
                    [content, comment] = replab.infra.splitComment(l(4:end));
                    content = strtrim(content);
                    comment = strtrim(comment);
                    if replab.compat.endsWith(content, '...')
                        tags(1, ind) = 'S';
                        content = content(1:end-3);
                    else
                        tags(1, ind) = 's';
                    end
                    lines{1, ind} = content;
                    comments{1, ind} = comment;
                    lineNumbers(1, ind) = i;
                else
                    [content, comment] = replab.infra.splitComment(l);
                    content = strtrim(content);
                    comment = strtrim(comment);
                    if replab.compat.endsWith(content, '...')
                        tags(1, ind) = 'O';
                        content = content(1:end-3);
                    else
                        tags(1, ind) = 'o';
                    end
                    lines{1, ind} = content;
                    comments{1, ind} = comment;
                    lineNumbers(1, ind) = i;
                end
                ind = ind + 1;
            end
            tags(1, ind) = '$';
            lines{1, ind} = '';
            comments{1, ind} = '';
            lineNumbers(1, ind) = 0;
            dtt = replab.infra.doctests.DocTestTokens(tags, lines, comments, lineNumbers);
        end

    end

end
