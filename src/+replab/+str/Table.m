classdef Table < replab.Str
% Helper to display tables on the REPL
%
% >>> replab.str.Table(cell(4,4), 2:3, 'llcc') % unreadable
%
% >>> replab.str.Table(cell(4,4), struct('omitRange', 2:3, 'colAlign', 'llcc'))
%
%
% >>> replab.str.Table(cell(4,4), 'omitRange', 2:3, 'colAlign', 'llcc')
%
% or
%
% >>> T = replab.str.Table(cell(4,4))
% >>> T.setColAlign(1:4, 'llcc')
% >>> T.setOmitRange(2:3)

% >>> T = replab.str.Table(elements);
% >>> T = T.withRowNames({'x1', 'x2', 'x3'}).withColumnNames(cols)

% why not? https://www.mathworks.com/matlabcentral/fileexchange/24093-cprintf-display-formatted-colored-text-in-the-command-wind
    properties
        rowSep % cell(1,nRows+1) of charstring: Row separators
        colSep % cell(1,nCols+1) of charstring: Column separators
               %
               %                                This includes a separator on the left of the first column and
               %                                on the right of the last column.
        colAlign % char(1,\*): Array of 'l' 'c' 'r'
        omitRange % integer(1,\*): Range of columns that can be omitted if space is lacking
                  %
                  %                This needs to be contiguous, and can be empty.
        elements % cell(\*,\*): Table contents
    end

    methods

        function self = Table(elements, varargin)
            s = struct(varargin{:});
            if isfield(s, 'omitRange')
                % ...
            end
            if isfield(s, 'colAlign')
                % ...
            end
        end

        function s = format(self, maxRows, maxColumns)
        % display code
        end

        function n = nColumns(self)
        end
        function n = nRows(self)
        end
        function setColSep(self, range, sep)
        % Sets the column separators for the columns inside the range
        %
        % Convention: the column to the left of the table is 0.
        %
        % Example:
        %   >>> T = replab.str.Table(elements);
        %   >>> T.setColSep(0:T.nColumns, '  ');
        % instead of replab.str.Table(elements, 'colSep', repmat({'  '}, 1, size(e

        end

        function setRowSep(self, range, sep)
        % same as setColSep

        % interesting problem
        %    |
        % ---?--- what to put in the '?'
        %    |
        end



        function addRowNames(self, rowNames)
        % Adds row names to the current table
        end

        function adColumnNames(self, colNames)
        % Adds column names to the current table
        end

    end

end

end
