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
        
        function flines = filteredLines(self)
        % Returns the documentation lines with some of the Sphinx formatting filtered out
            flines = cellfun(@(l) replab.infra.Doc.filterSphinxLine(l), self.lines, 'uniform', 0);
        end
        
        function dispFilteredLines(self, keyword, helpFunctionName, fullMode)
            flines = self.filteredLines;
            for i = 1:length(flines)
                if isempty(flines{i})
                    disp(' ');
                else
                    replab.infra.dispH(['   ', flines{i}], keyword, helpFunctionName, fullMode);
                end
            end
        end
        
        function b = isempty(self)
            b = true;
            for i = 1:length(self.lines)
                b = b && isempty(strtrim(self.lines{i}));
            end
        end

    end
    
    methods (Static)
        
        function l = filterSphinxLine(l)
        % Filters the Sphinx formatting out of a line
            
        % 1. we replace double backticks by single quotes
            l = strrep(l, '``', '''');
            % 2. we remove the plus '+' inside single backticks `bla`, and remove the backticks
            parts = strsplit(l, '`');
            for i = 2:2:length(parts)
                parts{i} = strrep(parts{i}, '+', '');
            end
            l = strjoin(parts, '');
        end
        
    end
    
end
