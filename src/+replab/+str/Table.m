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

    properties
        colAlign % char(1,\*): Array of 'l' 'c' 'r'
        omitRange % integer(1,\*): Range of columns that can be omitted if space is lacking, needs to be contiguous, can be empty
        elements % cell(\*,\*): Table contents
    end

    methods

        function self = Table(elements, varargin)
            s = struct(varargin{:});
            if isfield(s, 'omitRange')
                % ...
            end
            if isfield(s, 'colAlign')
                % ...x
            end
        end

        function s = format(self, maxRows, maxColumns)
        % display code
        end

        function newTable = withRowNames(self, rowNames)
        % Returns a new table with the added row names
        end

        function newTable = withColumnNames(self, colNames)
        % Returns a new table with the added column names
        end

    end

end
