classdef Table < replab.Str
%
% Displays tables nicely in MATLAB and Octave
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
               %                                This includes a separator above the first row and
               %                                below the last row.
        colSep % (cell(1,nColumns+1) of charstring): Column separators
               %
               %                                This includes a separator on the left of the first column and
               %                                on the right of the last column.
        colAlign % (char(1,nColumns)): Array of 'l' 'c' 'r' (default to 'c')
        omitRange % (integer(1,\*)): Range of columns that can be omitted if space is lacking
                  %
                  %                             This can be empty and must be a vector with no repeated values
        elements % (cell(nRows,nColumns)): Table contents (includes row and column names)
        columnLengths % (integer(1, nColumns)): Character width of each column in table
        rowName % (logical): Whether table has row names
        colName % (logical): Whether table has column names
    end

    
    
    methods

        function self = Table(elements, varargin)
        % 
        % Args:
        %   elements (cell(nRows, nColumns) or double(nRows, nColumns)): body of table
        %   varargin: list of 'name', value pairs to set table parameters
        %               'colAlign' column alignment character array
        %               'colSep' characters that separate columns
        %               'colName' cell array of column names
        %               'rowName' cell array of row names
        %               'omitRange' matrix of columns to omit first
        %               'rowSep' characters that separate rows
        %
        % Example:
        %   >>> T = replab.str.Table([11,2;100,4], 'colAlign', 'rlc', 'colSep', ' | ', ...
        %                 'colName', {'one', 'two'}, 'rowName', {1, 'second'}, ...
        %                 'omitRange', [1], 'rowSep', '-')
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
            self.elements = self.convertToChars(elements);
            
            elmtLens = max(cellfun(@length, self.elements), [], 1);
            sepLens = repmat(2, 1, self.nColumns + 1); % this is default
            self.columnLengths = [elmtLens, 0] + sepLens;
            
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
            
        end
        
        function disp(self)
            disp(self.format(replab.globals.strMaxRows, replab.globals.strMaxColumns))
        end

        function tbstr = format(self, maxRows, maxColumns)
        %
        % Formatting for the table display
        % Convention: - adds column separators and then row separators if given
        %             - if the table is wider than the display, omits first the columns
        %               given in omitRange and then works backwards from the
        %               second last column
        %             - if number of rows exceeds maxRows, omits rows backwards from
        %               second last row
        %
            dim = size(self.elements);
            char_arr = self.elements;
            
            % omit columns if table will be wider than maxColumns
            omitSymbol = ' ...';
            lens = self.columnLengths;
            colseps = self.colSep;
            len = sum(lens);
            spec = self.colAlign;
            if len > maxColumns
                ellipsisCol = cell(dim(1), 1);
                ellipsisCol{floor(dim(1)/2)} = omitSymbol;
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
                if len > maxColumns
                    colToOmit = dim(2) - 1;
                    len = len - lens(colToOmit) + length(omitSymbol);
                    while len > maxColumns 
                        colToOmit = colToOmit - 1;
                        len = len - lens(colToOmit);
                    end
                    char_arr(:, dim(2) - 1) = ellipsisCol;
                    char_arr(:, colToOmit:dim(2)-2) = [];
                    dim = size(char_arr);
                    colseps(dim(2) - 1) = {''};
                    colseps(colToOmit:dim(2)-2) = [];
                    spec(colToOmit:dim(2)-2) = [];
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
                    if spec(j) == 'l'
                        padding = elmt_lens(j) - length(elmt);
                        sized_elmt = [colseps{j}, elmt, repmat(' ', 1, padding)];
                    elseif spec(j) == 'c'
                        if rem(elmt_lens(j), 2) == 0
                            padding = floor((elmt_lens(j) - length(elmt)) / 2);
                        else
                            padding = ceil((elmt_lens(j) - length(elmt)) / 2);
                        end
                        sized_elmt = [colseps{j}, repmat(' ', 1, elmt_lens(j) - length(elmt) - padding), elmt, ...
                                      repmat(' ', 1, padding)];
                    elseif spec(j) == 'r'
                        padding = elmt_lens(j) - length(elmt);
                        sized_elmt = [colseps{j}, repmat(' ', 1, padding), elmt];
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
        end
        
        function n = nColumns(self)
        % 
        % number of columns (including row names)
            dim = size(self.elements);
            n = dim(2);
        end
        
        function n = nRows(self)
        %
        % number of rows (including column names)
            dim = size(self.elements);
            n = dim(1);
        end
        
        function setColSep(self, range, sep)
        %
        % Sets the column separators for the columns inside the range
        %
        % Convention: the column to the left of the table is 0.
        %
        % Example:
        %   >>> T = replab.str.Table({'one', 'two'});
        %   >>> T.setColSep(0:T.nColumns, ' | ');
        %   >>> T
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
            old_lengths = cellfun(@length, self.colSep);
            new_lengths = cellfun(@length, sep_array);
            self.columnLengths(range+1) = self.columnLengths(range+1) - old_lengths(range+1) + new_lengths;
            self.colSep(range+1) = sep_array;
        end

        function setRowSep(self, range, sep)
        % 
        % Sets the row separators for the rows inside the range
        %
        % Convention: the row at the top of the table is 0.
        %
        % Example:
        %   >>> T = replab.str.Table({'one'; 'two'});
        %   >>> T.setRowSep(0:T.nRows, '-');
        %   >>> T
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
        %
        % Set the alignment of each column in range
        % 'l' for left, 'c' for centre, and 'r' for right
        %
        % Example:
        %   >>> T = replab.str.Table({'one','two'; 'three','four'});
        %   >>> T.setAlign(1:T.nColumns, 'lr')
        %   >>> T
        %         one     two  
        %         three  four  
            % disregard characters after the range of the columns
            range(range > self.nColumns) = [];
            if length(aligns) > length(range)
                aligns = aligns(1:length(range));
            end
            self.colAlign(range) = aligns;
        end
        
        function setOmitRange(self, omitRange)
        %
        % Set columns that can be omitted if there is not enough space
        % 
        % Example: 
        %   >>> T = replab.str.Table({1:6,1:6,1:6,1:6,1:6,1:6;1,2,3,4,5,6});
        %   >>> T.setOmitRange([2:4])
        %   >>> disp(T.format(5, 75))
        %         [1, 2, 3, 4, 5, 6] ...  [1, 2, 3, 4, 5, 6]  [1, 2, 3, 4, 5, 6]  
        %                  1                       5                   6         
            self.omitRange = sort(omitRange);
        end

        function addRow(self, row, loc, sep)
        %
        % Adds row to the table
        %
        % Convention: loc = 0 means add a row above the table and 
        %             loc = nRows adds a row below the table
        %
        % Args:
        %   row (cell(1,nRows)): row to add to table
        %   loc (integer): location in table to add row
        %   sep (optional character array): row separator above new row
        % 
        % Example:
        %   >>> T = replab.str.Table([1,2,3]);
        %   >>> T.addRow({'one', 'two', 'three'}, T.nRows, '-')
        %   >>> T
        %          1    2     3  
        %       -------------------
        %         one  two  three 
            dim = size(self.elements);
            array = cell(dim(1) + 1, dim(2)); 
            array(loc + 1, :) = self.convertToChars(row);
            array(loc + 2:end, :) = self.elements(loc + 1:end, :);
            array(1:loc, :) = self.elements(1:loc, :);
            self.elements = array;
            if exist('sep', 'var')
                if ~isempty(self.rowSep)
                    self.rowSep = [self.rowSep(1:loc), {sep}, self.rowSep(loc+1:end)];
                else
                    self.rowSep = [cell(1, loc), {sep}, cell(1, self.nRows - loc)];
                end
            end
            elmtLens = max(cellfun(@length, self.elements));
            sepLens = cellfun(@length, self.colSep);
            self.columnLengths = [elmtLens, 0] + sepLens;
        end
        
        function addColumn(self, column, loc, sep, align)
        % adds row to the table
        %
        % Convention: loc = 0 means add a column to left of the table and 
        %             loc = T.nColumns adds a column to right of the table
        %
        % Args:
        %   column (cell(1,nColumns)): column to add to table
        %   loc (integer): location in table to add column
        %   sep (optional character array): column separator to left of new row
        % 
        % Example:
        %   >>> T = replab.str.Table([1;2;3]);
        %   >>> T.addColumn({'one', 'two', 'three'}, T.nColumns, ': ', 'l')
        %   >>> T
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
            colchar = self.convertToChars(column);
            array(:, loc + 1) = colchar;
            array(:, loc + 2:end) = self.elements(:, loc + 1:end);
            array(:, 1:loc) = self.elements(:, 1:loc);
            self.elements = array;
            self.colAlign = [self.colAlign(1:loc), align, self.colAlign(loc+1:end)];
            self.colSep = [self.colSep(1:loc), {sep}, self.colSep(loc+1:end)];
            columnwidth = max(cellfun(@length, colchar));
            self.columnLengths = [self.columnLengths(1:loc), ...
                                    columnwidth + length(sep), self.columnLengths(min(loc+1, end):end)];
        end

        function addColumnNames(self, colNames)
        %
        % Adds row names to the current table
        % Convention: should not include an entry above row names
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
        %
        % Adds column names to the current table
        % Convention: should not include an entry before column names
            if self.rowName
                self.elements = self.elements(:, 2:end);
                self.colAlign = self.colAlign(2:end);
                self.colSep = self.colSep(2:end);
                self.columnLengths = self.columnLengths(2:end);
            end
            if self.colName
                self.addColumn([{''}, rowNames], 0, '', 'c');
            else
                self.addColumn(rowNames, 0, '', 'c');
            end
            self.rowName = true;
        end
        
        function names = getColumnNames(self)
        % 
        % Returns:
        %   names (cell(1, nColumns-1)): column names if present or else
        %                                empty array
            if self.colName
                names = self.elements(1, 2:end);
            else
                names = {};
            end
        end
        
        function names = getRowNames(self)
        % 
        % Returns:
        %   names (cell(1, nRows-1)): row names if present or else
        %                                empty array
            if self.rowName
                names = self.elements(2:end, 1)';
            else
                names = {};
            end
        end


    end
    
    methods (Static)
        
        function chararray = convertToChars(elements)
            %
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
                    end
                    chararray{i, j} = charac;
                end
            end
        end
        
        function [ellipsisLocs, hiddenRange] = addEllipses(omitRange)
            %
            % puts ellipsis in the first of a range of omitted columns
            diff = omitRange(2:end) - omitRange(1:end-1);
            f = find(diff ~= 1);
            ellipsisLocs = omitRange([1, f + 1]);
            omitRange([1,f + 1]) = [];
            hiddenRange = omitRange;
        end
       
        
    end

end