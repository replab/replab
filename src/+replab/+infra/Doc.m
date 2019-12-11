classdef Doc < replab.Str
    
    properties
        lines % row cell array of charstring: Documentation lines
    end
    
    methods
        
        function self = Doc(lines)
            self.lines = lines;
        end
        
        function fl = firstLine(self)
            if ~isempty(self.lines)
                fl = self.lines{1};
            else
                fl = '';
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
        
        function dc = leftTrimmed(lines)
        % Constructs a DocComments instance with the left whitespace uniformly trimmed
        %
        % We count the amount of leading whitespace for each line, and computes the minimum
        % among all lines, then remove that number of whitespace characters from each line.
        %
        % We thus trim unnecessary leading whitespace while preserving the semantic information
        % the remaining whitespace contains.
            if ~isempty(lines)
                l = [];
                for i = 1:length(lines)
                    cl = lines{i};
                    if ~isempty(strtrim(cl))
                        token = regexp(cl, '^(\s*)', 'tokens', 'once');
                        if isempty(l)
                            l = length(token);
                        else
                            l = min(l, length(token));
                        end
                    end
                end
                if ~isempty(l)
                    lines = cellfun(@(x) x(l+1:end), lines, 'uniform', 0);
                end
            end
            dc = replab.infra.Doc(lines);
        end
        
    end
    
end
