classdef CodeTokens < replab.Str
% Stores the contents of a source code file
%
% We process the source code files line by line, i.e. each line is a token
% for the parser; whitespace is trimmed from the lines we read.
%
% We do not filter blank lines, so line numbers match (unlike DocTestParseState).
%
% Each line has an associated tag.
%
% Possible tags are:
%
% - ' ': The line contains only whitespace
% - '%': The line contains only a comment
% - '!': The line contains code
% - 'c: The line contains a 'c'lassdef statement
% - 'p': The line contains a 'p'roperties statement
% - 'm': The line contains a 'm'ethods statement
% - 'e': The line contains an 'e'numeration statement
% - 'f': The line contains a function statement
% - '>': The line contains an if/try/while/switch/parfor/spmd/for statement
% - '<: The line contains an end statement
% - '$': The end of file token
%
% See https://jayconrod.com/posts/65/how-to-build-a-parser-by-hand
% for another example of a hand-written parser
%
% However, DocTestParseState is closer to that example; here we don't allocate new parse states,
% we instead mutate the line position.

    properties
        filename % charstring: Filename or ``[]``
        lines % row cell array of charstring: Source code lines
        tags % charstring: Tag describing the line type, one char per tag
    end

    methods

        function self = CodeTokens(filename, lines, tags)
            self.filename = filename;
            self.tags = tags;
            self.lines = lines;
        end

        function id = sourceIdentifier(self)
        % Returns the source identifier, i.e. the name of the .m file without extension
            if isempty(self.filename)
                id = [];
            else
                [~,id,~] = fileparts(self.filename);
            end
        end

        function data = parse(self)
            switch self.peek(1)
              case 'c'
                data = replab.infra.ClassData.parse(self);
              case 'f'
                [pos data] = replab.infra.FunctionLikeData.parse(self, 1, []);
                if ~isequal(self.sourceIdentifier, data.name)
                    replab.infra.parseError(self, pos, 'Function declaration name %s does not match filename %s.m', ...
                                            data.name, self.sourceIdentifier);
                end
              otherwise
                replab.infra.parseError(self, 1, 'Invalid token at start of file %s.m', self.sourceIdentifier);
            end
        end

        function n = nLines(self)
        % Returns the number of actual lines, not counting the "end of file" added line
            n = length(self.lines) - 1;
        end

        function [nextPos tag line] = take(self, pos)
        % Consumes a line from the input, and returns the updated parser state
        %
        % ``tag`` and ``line`` correspond to the consumed line.
            tag = self.tags(pos);
            line = self.lines{pos};
            nextPos = pos + 1;
        end

        function [nextPos line] = expect(self, pos, expectedTag)
        % Equivalent to ``take``, conditioned on the consumed tag to be ``expectedTag``
            [nextPos tag line] = self.take(pos);
            if tag ~= expectedTag
                nextPos = [];
                line = [];
            end
        end

        function tag = peek(self, pos)
        % Peeks at the tag of the next line to be consumed, without consuming it
            tag = self.tags(pos);
        end

    end

    methods (Static)

        function ct = lex(filename, source)
        % Constructs a CodeTokens instance from source code lines
        %
        % Args:
        %   source (charstring): Source code
        %   lines (row cell array of charstring): Trimmed source code lines
        %
        % Returns:
        %   :class:`+replab.+infra.CodeTokens`: A fresh CodeTokens instance
            lines = cellfun(@strtrim, strsplit(source, '\n', 'CollapseDelimiters', false), 'uniform', 0);
            n = length(lines);
            tags = blanks(n+1);
            tags(n+1) = '$';
            lines{n+1} = '';
            for i = 1:n
                line = lines{i};
                if isempty(line)
                    tags(i) = ' ';
                elseif line(1) == '%'
                    tags(i) = '%';
                else
                    word = regexp(line, '^\w+', 'match');
                    if isempty(word)
                        tags(i) = '!';
                    else
                        switch word{1}
                          case 'classdef'
                            tags(i) = 'c';
                          case 'properties'
                            tags(i) = 'p';
                          case 'methods'
                            tags(i) = 'm';
                          case 'enumeration'
                            tags(i) = 'e';
                          case 'function'
                            tags(i) = 'f';
                          case {'if', 'try', 'while', 'switch', 'parfor', 'spmd', 'for'}
                            tags(i) = '>';
                          case 'end'
                            tags(i) = '<';
                          otherwise
                            tags(i) = '!';
                        end
                    end
                end
            end
            ct = replab.infra.CodeTokens(filename, lines, tags);
        end

        function ct = fromSource(source)
        % Constructs a CodeTokens instance from the source code
        %
        % Args:
        %   source (charstring): Source code a single char string
        %
        % Returns:
        %   :class:`+replab.+infra.CodeTokens`: A fresh CodeTokens instance
            ct = replab.infra.CodeTokens.lex([], contents);
        end

        function ct = fromFile(filename)
        % Constructs a CodeTokens instance from the contents of the given filename
        %
        % Args:
        %   filename (charstring): Name of the source file to read and lex
        %
        % Returns:
        %   :class:`+replab.+infra.CodeTokens`: A fresh CodeTokens instance
            contents = fileread(filename);
            ct = replab.infra.CodeTokens.lex(filename, contents);
        end

    end

end
