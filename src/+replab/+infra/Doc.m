classdef Doc < replab.Str
% Contents of a documentation block comment

    properties
        parent % `.SourceElement`: Member containing this documentation
        lines % row cell array of charstring: Trimmed documentation lines
        lineNumbers % row integer vector: Line numbers
    end

    methods

        function self = Doc(parent, rawLines, lineNumbers)
            self.parent = parent;
            self.lines = replab.infra.uniformLeftTrim(rawLines);
            self.lineNumbers = lineNumbers;
        end

        function c = content(self)
        % Returns the lines joined with newlines
            c = strjoin(self.lines, char(10));
        end

        function n = nLines(self)
        % Returns the number of lines in this documentation comment block
            n = length(self.lines);
        end

        function l = line(self, i)
        % Returns the i-th line of this documentation comment block
            l = self.lines{i};
        end

        function fl = firstLine(self)
        % Returns the first line of the documentation comment block
        %
        % Care is taken to return an empty string if the comment block is empty,
        % and to add trailing ``...`` if for some reason the first paragraph
        % spans multiple lines.
            switch self.nLines
              case 0
                fl = '';
              case 1
                fl = self.line(1);
              otherwise
                if isempty(self.line(2))
                    fl = self.line(1);
                else
                    fl = [self.line(1) ' ...'];
                end
            end
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Str(self);
            if ~isempty(lines)
                names{1,end+1} = 'firstLine';
                values{1,end+1} = self.firstLine;
            end
        end

        function b = isempty(self)
            b = true;
            for i = 1:length(self.lines)
                b = b && isempty(strtrim(self.lines{i}));
            end
        end

    end

end
