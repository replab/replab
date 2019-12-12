classdef DocTestParseState < replab.Str
% Stores the contents and the parsing position (=state) of a doctest block 
%
% We process the doctest block line by line, i.e. each line is a token
% for the parser; whitespace is trimmed from the lines we read.
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
        lines % row cell array of charstring: Content lines
        comments % row cell array of charstring: Comments
        tags % row cell array of charstring: Tag describing the line type
        pos % integer: Current line position
    end
    
    methods
        
        function self = DocTestParseState(lines, comments, tags, pos)
            self.tags = tags;
            self.lines = lines;
            self.comments = comments;
            self.pos = pos;
        end
        
        function [nextParseState tag line comment] = take(self)
        % Consumes a line from the input, and returns the updated parser state
        %
        % `tag` and `line` correspond to the consumed line.
            tag = self.tags{self.pos};
            if isequal(tag, 'EOF')
                line = [];
                comment = [];
            else
                line = self.lines{self.pos};
                comment = self.comments{self.pos};
            end
            nextParseState = replab.infra.DocTestParseState(self.lines, self.comments, self.tags, self.pos + 1);
        end
        
        function [nextParseState line comment] = expect(self, expectedTag)
        % Equivalent to `take`, conditioned on the consumed tag to be `expectedTag`
            [nextParseState tag line comment] = self.take;
            if ~isequal(tag, expectedTag)
                nextParseState = [];
                line = [];
                comment = [];
            end
        end
        
    end
    
    methods (Static)

        function [code comment] = splitComment(l)
        % Splits a line around the '%' character that denotes a comment, ignoring the '%' inside single quoted strings
        %
        % Works also when the source line contains no such character.
        %
        % Args:
        %   l (charstring): Source line to split
        %
        % Returns
        % -------
        %   code:
        %     charstring: The code line before the '%'
        %   comment:
        %     charstring: The comment after the '%' (if there is no comment, this return value is empty)
            q = (l == ''''); % single quotes
            insideComment = logical(bitand(cumsum(q), 1)); % boolean mask of comments, excluding final single quotes
            commentStart = and(~insideComment, l == '%');
            pos = find(commentStart);
            if ~isempty(pos)
                code = l(1:pos(1)-1);
                comment = l(pos(1)+1:end);
            else
                code = l;
                comment = '';
            end
        end

        function ps = fromDocTestBlock(inputLines)
        % Constructs a ParseState instance from a doctest block lines
        %
        % Args:
        %   lines (row cell array of charstring): Doctest block lines
        %
        % Returns:
        %   :class:`replab.infra.ParseState`: A fresh ParseState instance
            lines = {};
            comments = {};
            tags = {};
            for i = 1:length(inputLines)
                l = strtrim(inputLines{i});
                if isempty(l) || l(1) == '%'
                    % omit that line
                elseif length(l) >= 3 && isequal(l(1:3), '>>>')
                    tags{1,end+1} = 'START';
                    [code comment] = replab.infra.DocTestParseState.splitComment(l(4:end));
                    lines{1,end+1} = strtrim(code);
                    comments{1,end+1} = strtrim(comment);
                elseif length(l) >= 3 && isequal(l(1:3), '...')
                    tags{1,end+1} = 'CONT';
                    [code comment] = replab.infra.DocTestParseState.splitComment(l(4:end));
                    lines{1,end+1} = strtrim(code);
                    comments{1,end+1} = strtrim(comment);
                else
                    tags{1,end+1} = 'OUT';
                    lines{1,end+1} = l;
                    comments{1,end+1} = '';
                end
            end
            tags{1,end+1} = 'EOF';
            ps = replab.infra.DocTestParseState(lines, comments, tags, 1);
        end
        
    end
    
end
