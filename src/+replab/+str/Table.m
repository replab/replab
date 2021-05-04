classdef Table < replab.Str
% Displays tables nicely in MATLAB and Octave
%
% Can take in any cell array or matrix
%
%   >>> replab.str.Table([1,20;300,2])
%         1   20
%        300   2
%   >>> replab.str.Table({'hello', 'world'; 'world', 'hi'})
%         hello  world
%         world   hi

    properties
        rowSep % (cell(1,nRows+1) of charstring): Row separators
               %
               %                                  This includes a separator above the first row and
               %                                  below the last row.
        colSep % (cell(1,nColumns+1) of charstring): Column separators
               %
               %                                  This includes a separator on the left of the first column and
               %                                  on the right of the last column.
        colAlign % (char(1,nColumns)): Array of 'l' 'c' 'r' (default to 'c')
        omitRange % (integer(1,\*)): Range of columns that can be omitted if space is lacking
                  %
                  %                  This can be empty and must be a vector with no repeated values
        elements % (cell(nRows,nColumns)): Table contents (includes row and column names)
        rowName % (logical): Whether table has row names
        colName % (logical): Whether table has column names
        title % (charstring): String to show above table (shows nothing if empty)
    end

    methods

        function self = Table(elements, varargin)
        % Creates a Table object
        %
        % The constructor can take a variable number of arguments after ``elements``, which are name-value pairs.
        %
        % - 'colAlign' column alignment charstring. Must be either one character that is repeated
        %              for all columns or charstring of length nColumns (including column of row names if given)
        % - 'rowSep' characters that separate rows. Charstring given will repeated to fill the row
        % - 'colSep' characters that separate columns. This will be used for all columns but individual
        %            columns can be set with .setColSep
        % - 'rowName' cell(1, nRows) array of row names
        % - 'colName' cell(1, nColumns) array of column names
        % - 'omitRange' integer(1,\*) vector of columns to omit first
        % - 'title' charstring of table title
        %
        % Args:
        %   elements (cell(nRows, nColumns) or double(nRows, nColumns)): body of table
        %   varargin: list of 'name', value pairs to set table parameters
        %
        % Example:
        %   >>> T = replab.str.Table([11,2;100,4], 'colAlign', 'rlc', 'colSep', ' | ', ...
        %                 'colName', {'one', 'two'}, 'rowName', {1, 'second'}, ...
        %                 'omitRange', [1], 'rowSep', '-', 'title', 'Table 1')
        %       Table 1
        %       ------------------------
        %        |        | one | two |
        %       ------------------------
        %        |      1 | 11  |  2  |
        %       ------------------------
        %        | second | 100 |  4  |
        %       ------------------------
        %
            if ~iscell(elements) && ismatrix(elements)
                elements = num2cell(elements);
            elseif ~iscell(elements)
                err('Error: elements must be matrix or cell array')
            end
            dim = size(elements);
            self.elements = elements;

            self.colName = false;
            self.rowName = false;

            % use variables if given and otherwise default settings
            s = struct(varargin{:});
            % variables with default settings
            if isfield(s, 'colAlign')
                if length(s(1).colAlign) == self.nColumns + 1 && isfield(s, 'colName')
                    self.setAlign(1:self.nColumns, s(1).colAlign(2:end))
                else
                    self.setAlign(1:self.nColumns, s(1).colAlign)
                end
            else
                self.setAlign(1:self.nColumns, repmat('c', 1, dim(2)))
            end
            if isfield(s, 'colSep')
                self.setColSep(0:self.nColumns, s(1).colSep);
            else
                self.setColSep(0:self.nColumns, '  ');
            end

            % optional variables to set
            if isfield(s, 'colName')
                self.addColumnNames(cellfun(@(x) s(x).colName, num2cell(1:dim(2)), ...
                                                'UniformOutput', false))
                self.colName = true;
            end
            if isfield(s, 'rowName')
                self.addRowNames(cellfun(@(x) s(x).rowName, num2cell(1:dim(1)), ...
                                          'UniformOutput', false))
                self.rowName = true;
                if isfield(s, 'colAlign') && length(s(1).colAlign) == self.nColumns
                    self.setAlign(1, s(1).colAlign(1))
                end
                if isfield(s, 'colSep')
                    self.setColSep(0, s(1).colSep)
                end
            end
            if isfield(s, 'omitRange')
                self.setOmitRange(s(1).omitRange);
            end
            if isfield(s, 'rowSep')
                self.setRowSep(0:self.nRows, s(1).rowSep);
            end
            if isfield(s, 'title')
                self.addTitle(s(1).title)
            end

        end

        function disp(self)
            disp(self.format(replab.globals.strMaxRows, replab.globals.strMaxColumns))
        end

        function s = headerStr(self)
            dim = size(self.elements);
            if self.colName
                dim(1) = dim(1) - 1;
            end
            if self.rowName
                dim(2) = dim(2) - 1;
            end
            s = sprintf([num2str(dim(1)), ' x ', num2str(dim(2)), ' Table']);
        end

        function [tbstr, truncated] = format(self, maxRows, maxColumns)
        % Formatting for the table display
        %
        % Convention: - use maxRows = {[],Inf} or maxColumns = {[],Inf} to have no restriction
        %               on output size
        %             - adds column separators and then row separators if given
        %             - if the table is wider than the display, omits first the columns
        %               given in omitRange and then works backwards from the
        %               second last column
        %             - if number of rows exceeds maxRows, omits rows backwards from
        %               second last row
        %
        % Args:
        %   maxRows (integer): maximum number of rows (including rows of separators)
        %   maxColumns (integer): maximum number of columns (where column
        %   refers to character on the display screen)
        %
        % Returns:
        %   tbstr (charstring): string with ``\n`` separators for display
        %   truncated (logical): whethere part of the table was omitted
            dim = size(self.elements);
            if nargin < 3 || isempty(maxColumns)
                maxColumns = Inf;
            end
            if nargin < 2 || isempty(maxRows)
                maxRows = Inf;
            end
            % to replace align, make sure that two column tables with
            % strings in the first column won't have a column omitted
            if dim(2) == 2 && all(cellfun(@ischar, self.elements(:, 1)))
                nameLen = max(cellfun(@length, self.elements(:,1)));
                char_arr = self.convertToChars(self.elements, maxColumns-nameLen-2);
            else
                char_arr = self.convertToChars(self.elements, maxColumns);
            end


            % omit columns if table will be wider than maxColumns
            omitSymbol = ' ...';
            truncated = false;
            elmtLens = max(cellfun(@length, char_arr), [], 1);
            sepLens = cellfun(@length, self.colSep);
            lens = [elmtLens,0] + sepLens;
            colseps = self.colSep;
            len = sum(lens);
            spec = self.colAlign;
            if len > maxColumns
                truncated = true;
                ellipsisCol = cell(dim(1), 1);
                ellipsisCol{max(floor(dim(1)/2), 1)} = omitSymbol;
                if ~isempty(self.omitRange)
                    [dots, hide] = self.addEllipses(self.omitRange);
                    char_arr(:, dots) = repmat(ellipsisCol, 1, length(dots));
                    char_arr(:, hide) = [];
                    dim = size(char_arr);
                    lens(dots) = 3;
                    lens(hide) = [];
                    colseps(dots) = {''};
                    colseps(hide) = [];
                    spec(hide) = [];
                    len = sum(lens);
                end
                % start omitting columns backwards from second last column
                if len > maxColumns && dim(2) > 1
                    colToOmit = dim(2) - 1;
                    len = len - lens(colToOmit) + length(omitSymbol);
                    while len > maxColumns && colToOmit > 2
                        colToOmit = colToOmit - 1;
                        len = len - lens(colToOmit);
                    end
                    char_arr(:, dim(2) - 1) = ellipsisCol;
                    if dim(2) > 2
                        char_arr(:, colToOmit:dim(2)-2) = [];
                    end
                    dim = size(char_arr);
                    colseps(dim(2) - 1) = {''};
                    if dim(2) > 2
                        colseps(colToOmit:dim(2)-2) = [];
                        spec(colToOmit:dim(2)-2) = [];
                    end
                end
            end

            % omit rows if number of rows is bigger than maxRows
            if ~isempty(self.rowSep)
                rowseps = self.rowSep;
            end
            emptyRows = [];
            if ~isempty(self.rowSep)
                emptyRows = cellfun(@isempty, rowseps);
            end
            omitRows = false;
            if dim(1) + length(self.rowSep) - sum(emptyRows) > maxRows
                truncated = true;
                if ~isempty(self.rowSep)
                    nrows = 2 + ~emptyRows(end);
                    i = 1;
                    while nrows < maxRows
                        nrows = nrows + 1 + ~emptyRows(i);
                        i = i + 1;
                    end
                    char_arr(i:end-1, :) = [];
                    rowseps(i:end-2) = [];
                    omitRows = true;
                else
                    nremove = dim(1) - maxRows;
                    char_arr(end-nremove:end-1, :) = [];
                    omitRows = true;
                end
            end

            % add padding to the characters to make all elements same length
            dim = size(char_arr);
            elmt_lens = max(cellfun(@length, char_arr), [], 1);
            for i = 1:dim(1)
                for j = 1:dim(2)
                    elmt = char_arr{i, j};
                    csj = colseps{j};
                    if isempty(csj)
                        csj = '';
                    end
                    if spec(j) == 'l'
                        padding = elmt_lens(j) - length(elmt);
                        sized_elmt = [csj, elmt, repmat(' ', 1, padding)];
                    elseif spec(j) == 'c'
                        if rem(elmt_lens(j), 2) == 0
                            padding = floor((elmt_lens(j) - length(elmt)) / 2);
                        else
                            padding = ceil((elmt_lens(j) - length(elmt)) / 2);
                        end
                        sized_elmt = [csj, repmat(' ', 1, elmt_lens(j) - length(elmt) - padding), elmt, ...
                                      repmat(' ', 1, padding)];
                    elseif spec(j) == 'r'
                        padding = elmt_lens(j) - length(elmt);
                        sized_elmt = [csj, repmat(' ', 1, padding), elmt];
                    end
                    char_arr{i, j} = sized_elmt;

                end
            end

            if ~isempty(colseps{end})
                char_arr(:, end+1) = mat2cell(repmat(colseps{end}, dim(1), 1), ones(1, dim(1)));
            end
            if ~isempty(self.rowSep)
                nrows = dim(1) + length(rowseps);
                tbarray = cell(nrows, 1);
                for i = 1:nrows
                    if rem(i, 2) == 0
                        tbarray{i} = strjoin(char_arr(i/2, :), '');
                    else
                        rowchar = rowseps{(i+1)/2};
                        if ~isempty(rowchar)
                            new_row = repmat(rowchar, 1, floor(sum(lens)/length(rowchar)));
                            if length(new_row) ~= sum(lens)
                                tbarray{i} = [new_row, rowchar(sum(lens) - length(new_row))];
                            else
                                tbarray{i} = new_row;
                            end
                        end
                    end
                end
                if omitRows
                      tbarray(end - 2) = {'...'};
                end
                tbarray(cellfun(@isempty, tbarray)) = [];
            else
                tbarray = cell(dim(1), 1);
                for i = 1:dim(1)
                    tbarray{i} = strjoin(char_arr(i, :), '');
                end
                if omitRows
                    tbarray(end - 1) = {'...'};
                end
            end
            tbstr = strjoin(tbarray, '\n');
            if ~isempty(self.title)
                tbstr = strjoin({self.title, tbstr}, '\n');
            end
        end

        function n = nColumns(self)
        % Returns the number of columns (including row names)
        %
        % Returns:
        %   n (integer): Number of columns
            dim = size(self.elements);
            n = dim(2);
        end

        function n = nRows(self)
        % Returns the number of rows (including column names)
        %
        % Returns:
        %   n (integer): Number of columns
            dim = size(self.elements);
            n = dim(1);
        end

        function setColSep(self, range, sep)
        % Sets the column separators for the columns inside the range
        %
        % Convention: the column to the left of the table is 0.
        %
        % Args:
        %   range (integer(\*)): vector of column separator locations
        %   sep (charstring): characters that separate columns in range
        %
        % Example:
        %   >>> T = replab.str.Table({'one', 'two'});
        %   >>> T.setColSep(0:T.nColumns, ' | ');
        %   >>> T
        %       T =
        %        | one | two |
        %
            if isempty(self.colSep)
                self.colSep = mat2cell(repmat('  ', self.nColumns+1, 1), ones(self.nColumns+1,1))';
            end
            % Don't let the number of separators exceed the table width
            range(range > self.nColumns) = [];
            if isempty(sep)
                sep_array = cell(1, length(range));
            else
                sep_array = mat2cell(repmat(sep, length(range), 1), ones(length(range),1))';
            end
            self.colSep(range+1) = sep_array;
        end

        function setRowSep(self, range, sep)
        % Sets the row separators for the rows inside the range
        %
        % Convention: the row at the top of the table is 0.
        %
        % Args:
        %   range (integer(\*)): vector of row separator locations
        %   sep (charstring): characters to be repeated to separate rows in range
        %
        % Example:
        %   >>> T = replab.str.Table({'one'; 'two'});
        %   >>> T.setRowSep(0:T.nRows, '-');
        %   >>> T
        %      T =
        %       -------
        %         one
        %       -------
        %         two
        %       -------
        %
            if isempty(self.rowSep)
                self.rowSep = cell(1, self.nRows + 1);
            end
            % Don't let the number of separators exceed the table width
            range(range > self.nRows) = [];
            if isempty(sep)
                self.rowSep(range+1) = cell(1, length(range));
            else
                self.rowSep(range+1) = mat2cell(repmat(sep, length(range), 1), ones(length(range),1));
            end
        end

        function setAlign(self, range, aligns)
        % Set the alignment of each column in range
        %
        % Args:
        %   range (integer(1,n)): vector of column positions in the same order as aligns
        %   aligns (char(1,n)):  'l' for left, 'c' for centre, and 'r' for right
        %
        % Example:
        %   >>> T = replab.str.Table({'one','two'; 'three','four'});
        %   >>> T.setAlign(1:T.nColumns, 'lr');
        %   >>> T
        %       T =
        %         one     two
        %         three  four
        %
            range(range > self.nColumns) = []; % disregard characters after the range of the columns
            if length(aligns) > length(range)
                aligns = aligns(1:length(range));
            elseif length(aligns) < length(range)
                repalign = repmat(aligns, 1, floor(length(range)/length(aligns)));
                aligns = [repalign, aligns(1:length(range)-length(repalign))];
            end
            self.colAlign(range) = aligns;
        end

        function setOmitRange(self, omitRange)
        % Set columns that can be omitted if there is not enough space
        %
        % Args:
        %   omitRange (integer(\*)): range of columns that will be omitted first if
        %                            columns don't fit
        %
        % Example:
        %   >>> T = replab.str.Table({1:6,1:6,1:6,1:6,1:6,1:6;1,2,3,4,5,6});
        %   >>> T.setOmitRange([2:4]);
        %   >>> disp(T.format(5, 75))
        %         [1, 2, 3, 4, 5, 6] ...  [1, 2, 3, 4, 5, 6]  [1, 2, 3, 4, 5, 6]
        %                  1                       5                   6
        %
            self.omitRange = sort(omitRange);
        end

        function addRow(self, row, loc, sep)
        % Adds row to the table
        %
        % Convention: loc = 0 means add a row above the table and
        %             loc = nRows adds a row below the table
        %
        % Args:
        %   row (cell(1,nRows)): Row to add to table
        %   loc (integer): Location in table to add row
        %   sep (charstring, optional): Row separator above new row,
        %                               default to empty ''
        %
        % Example:
        %   >>> T = replab.str.Table([1,2,3]);
        %   >>> T.addRow({'one', 'two', 'three'}, T.nRows, '-');
        %   >>> T
        %       T =
        %          1    2     3
        %       -------------------
        %         one  two  three
            dim = size(self.elements);
            array = cell(dim(1) + 1, dim(2));
            array(loc + 1, :) = row;
            array(loc + 2:end, :) = self.elements(loc + 1:end, :);
            array(1:loc, :) = self.elements(1:loc, :);
            self.elements = array;
            if ~exist('sep', 'var')
                sep = '';
            end
            if ~isempty(self.rowSep)
                self.rowSep = [self.rowSep(1:loc), {sep}, self.rowSep(loc+1:end)];
            else
                self.rowSep = [cell(1, loc), {sep}, cell(1, self.nRows - loc)];
            end
        end

        function addColumn(self, column, loc, sep, align)
        % adds column to the table
        %
        % Convention: loc = 0 means add a column to left of the table and
        %             loc = T.nColumns adds a column to right of the table
        %
        % Args:
        %   column (cell(1,nColumns)): Column to add to table
        %   loc (integer): Location in table to add column
        %   sep (charstring, optional): Column separator to left of new row
        %
        % Example:
        %   >>> T = replab.str.Table([1;2;3]);
        %   >>> T.addColumn({'one', 'two', 'three'}, T.nColumns, ': ', 'l');
        %   >>> T
        %       T =
        %         1: one
        %         2: two
        %         3: three
            if ~exist('sep', 'var')
                sep = '  ';
            end
            if ~exist('align', 'var')
                align = 'c';
            end
            dim = size(self.elements);
            array = cell(dim(1), dim(2) + 1);
            array(:, loc + 1) = column;
            array(:, loc + 2:end) = self.elements(:, loc + 1:end);
            array(:, 1:loc) = self.elements(:, 1:loc);
            self.elements = array;
            self.colAlign = [self.colAlign(1:loc), align, self.colAlign(loc+1:end)];
            self.colSep = [self.colSep(1:loc), {sep}, self.colSep(loc+1:end)];
        end

        function addColumnNames(self, colNames)
        % Adds column names to the current table
        %
        % Convention: should not include an entry above row names
        %
        % Args:
        %   colNames (cell(1,nColumns)): cell array of column names
            if self.colName
                self.elements = self.elements(2:end,:);
                self.rowSep = self.rowSep(2:end);
            end
            if self.rowName
                self.addRow([{''}, colNames], 0, '');
            else
                self.addRow(colNames, 0, '');
            end
            self.colName = true;
        end

        function addRowNames(self, rowNames)
        % Adds row names to the current table
        %
        % Convention: should not include an entry before column names
        %
        % Args:
        %   rowNames (cell(1,nRows)): cell array of row names
            if self.rowName
                self.elements = self.elements(:, 2:end);
                self.colAlign = self.colAlign(2:end);
                self.colSep = self.colSep(2:end);
            end
            if self.colName
                self.addColumn([{''}, rowNames], 0, '', 'c');
            else
                self.addColumn(rowNames, 0, '', 'c');
            end
            self.rowName = true;
        end

        function addTitle(self, title)
        % Adds title to the table
        %
        % Args:
        %   title (charstring): charstring to display above the table
            if ischar(title)
                self.title = title;
            end
        end

        function names = getColumnNames(self)
        % Returns the column names
        %
        % Returns:
        %   names (cell(1, nColumns-1)): column names if present or else empty array
            if self.colName
                if self.rowName
                    names = self.elements(1, 2:end);
                else
                    names = self.elements(1, 1:end);
                end
            else
                names = {};
            end
        end

        function names = getRowNames(self)
        % Returns the row names
        %
        % Returns:
        %   names (cell(1, nRows-1)): row names if present or else empty array
            if self.rowName
                if self.rowName
                    names = self.elements(2:end, 1)';
                else
                    names = self.elements(1:end, 1)';
                end
            else
                names = {};
            end
        end

        function row = row(self, loc)
        % Returns the row at the given location
        %
        % Convention: - returns as a vector only if all entries are numeric
        %             - loc does not include the row of column names
        %             - row will not include the row name
        %
        % Args:
        %   loc (integer): row number
        %
        % Returns:
        %   row ({cell(1,\*), double(1,\*)}): cell array of row entries or vector of row
        %                                     entries if all are numeric
            if self.colName
                loc = loc + 1;
            end
            if self.rowName
                row = self.elements(loc, 2:end);
            else
                row = self.elements(loc, 1:end);
            end
            if all(cellfun(@isnumeric, row))
                row = cell2mat(row);
            end
        end

        function col = column(self, loc)
        % Returns the column at the given location
        %
        % Convention: - returns as a vector only if all entries are numeric
        %             - loc does not include the column of row names
        %             - col will not include the column name
        %
        % Args:
        %   loc (integer): row number
        %
        % Returns:
        %   col ({cell(1,\*), double(1,\*)}): cell array of column entries or vector of column
        %                                     entries if all are numeric
            if self.rowName
                loc = loc + 1;
            end
            if self.colName
                col = self.elements(2:end, loc);
            else
                col = self.elements(1:end, loc);
            end
            if all(cellfun(@isnumeric, col))
                col = cell2mat(col);
            end
        end


    end

    methods (Static)

        function chararray = convertToChars(elements, maxColumns)
            % converts all of elements to characters
            %
            % Args:
            %    elements (cell array or matrix)
            %
            % Returns:
            %    chararray (cell array of strings)

            if ~iscell(elements) && ismatrix(elements)
                elements = num2cell(elements);
            end
            dim = size(elements);
            chararray = cell(dim);
            for i = 1:dim(1)
                for j = 1:dim(2)
                    elmt = elements{i, j};
                    if ischar(elmt)
                        charac = elmt;
                    elseif isstring(elmt) % will only be true in MATLAB
                        charac = convertStringsToChars(elmt);
                    else
                        charac = replab.shortStr(elmt);
                        if length(charac) > maxColumns && ismatrix(elmt)
                            charac = replab.headerStr(elmt);
                        end
                    end
                    if length(charac) > maxColumns
                        chararray{i, j} = [charac(1:maxColumns-3), '...'];
                    else
                        chararray{i, j} = charac;
                    end
                end
            end
        end

        function [ellipsisLocs, hiddenRange] = addEllipses(omitRange)
            % puts ellipsis in the first of a range of omitted columns
            diff = omitRange(2:end) - omitRange(1:end-1);
            f = find(diff ~= 1);
            ellipsisLocs = omitRange([1, f + 1]);
            omitRange([1,f + 1]) = [];
            hiddenRange = omitRange;
        end

    end

end