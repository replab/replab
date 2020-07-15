classdef Table < replab.Str
% Helper to display tables on the REPL
%
% >>> replab.str.Table(cell(4,4), 2:3, 'llcc') % unreadable

    properties
        rowSep % cell(1,nRows+1) of charstring: Row separators
        colSep % cell(1,nCols+1) of charstring: Column separators
               %
               %                                This includes a separator on the left of the first column and
               %                                on the right of the last column.
        colAlign % char(1,nCols): Array of 'l' 'c' 'r' (default to 'c')
        omitRange % integer(1,\*): Range of columns that can be omitted if space is lacking
                  %
                  %                This can be empty.
        elements % cell(\*,\*): Table contents
        columnLengths % integer(1, nCols): character width of each column in table
        rowName % {1, 0}: whether table has row names
        colName % {1, 0}: whether table has column names
    end

    
    
    methods

        function self = Table(elements, varargin)
            if ~iscell(elements) && ismatrix(elements)
                elements = num2cell(elements);
            elseif ~iscell(elements)
                err('Error: elements must be matrix or cell array')
            end
            self.colName = 0;
            self.rowName = 0;
            s = struct(varargin{:});
            if isfield(s, 'Cols')
                self.addColumnNames(s.Cols)
                self.colName = 1;
            end
            if isfield(s, 'Rows')
                self.addRowNames(s.Rows)
                self.rowName = 1;
            end
            dim = size(elements);
            self.elements = self.convertToChars(elements);
            
            % use variables if given and otherwise default settings  
            if isfield(s, 'omitRange')
                self.omitRange = s.omitRange;
            end
            if isfield(s, 'colAlign')
                self.colAlign = s.colAlign;
            else
                self.colAlign = repmat('c', 1, dim(2));
            end
            if isfield(s, 'rowSep')
                self.rowSep = s.rowSep;
            end
            if isfield(s, 'colSep')
                self.colSep = s.colSep;
            else
                self.colSep = mat2cell(repmat('  ', dim(2)+1, 1), ones(1, dim(2)+1))';
            end
            
            % determine width of each column
            elmtLens = max(cellfun(@length, elements));
            sepLens = cellfun(@length, self.colSep);
            self.columnLengths = [elmtLens, 0] + sepLens;
        end
        
        function disp(self)
            disp(self.format(replab.globals.strMaxRows, replab.globals.strMaxColumns))
        end

        function tbstr = format(self, maxRows, maxColumns)
            dim = size(self.elements);
            char_arr = self.elements;
            
            % omit columns if table will be wider than maxColumns
            lens = self.columnLengths;
            colseps = self.colSep;
            len = sum(lens);
            spec = self.colAlign;
            if len > maxColumns
                ellipsisCol = cell(1, dim(2));
                ellipsisCol{floor(dim(2)/2)} = '...';
                if ~isempty(self.omitRange)
                    [dots, hide] = self.addEllipses(self.omitRange);
                    char_arr(:, dots) = ellipsisCol;
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
                    len = len - lens(colToOmit) + 3;
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
                nrows = dim(1);
            else
                nrows = dim(1) + length(self.rowSep);
            end
            if nrows > maxRows
                if ~isempty(self.rowSep)
                    nremove = ceil((maxRows - nrows - 1)/2);
                    rowseps(end-nremove:end-1) = [];
                    rowseps(end-1) = {''};
                    char_arr(end-nremove:end-1, :) = [];
                    char_arr(end-1) = {'...', cell(1, dim(2))};
                else
                    nremove = maxRows - nrows - 2;
                    char_arr(end-nremove:end-1, :) = [];
                    char_arr(end-nremove-1) = {'...', cell(1, dim(2))};
                end
            end
            
            % add padding to the characters to make all elements same length 
            dim = size(char_arr);
            elmt_lens = max(cellfun(@length, char_arr));
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
                nrows = dim(1) + length(self.rowSep);
                tbarray = cell(nrows, 1);
                for i = 1:nrows
                    if rem(i, 2) == 0
                        disp(char_arr(i/2, :))
                        tbarray{i} = strjoin(char_arr(i/2, :), '');
                    else
                        rowchar = rowseps{i/2-1};
                        new_row = repmat(rowchar, 1, floor(sum(lens)/length(rowchar)));
                        tbarray{i} = [new_row, rowchar(sum(lens) - length(new_row))];
                    end
                end
            else
                tbarray = cell(dim(1), 1);
                testarray = char_arr(i, :);
                for i = 1:dim(1)
                    tbarray{i} = strjoin(char_arr(i, :), '');
                end
            end
            tbstr = strjoin(tbarray, '\n');
        end
        
        function n = nColumns(self)
            dim = size(self.elements);
            n = dim(2);
        end
        
        function n = nRows(self)
            dim = size(self.elements);
            n = dim(1);
        end
        
        function setColSep(self, range, sep)
        % Sets the column separators for the columns inside the range
        %
        % Convention: the column to the left of the table is 0.
        %
        % Example:
        %   >>> T = replab.str.Table(elements);
        %   >>> T.setColSep(0:T.nColumns, '  ');
        % 
            if isempty(sep)
                self.colSep(range+1) = cell(length(range), 1);
            else
                self.colSep(range+1) = mat2cell(repmat(sep, length(range), 1), ones(length(range),1));
            end
        end

        function setRowSep(self, range, sep)
        % Sets the row separators for the rows inside the range
        %
        % Convention: the row at the top of the table is 0.
        %
        % Example:
        %   >>> T = replab.str.Table(elements);
        %   >>> T.setRowSep(0:T.nRows, '  ');
        % 
            self.rowSep(range+1) = mat2cell(repmat(sep, length(range), 1), ones(length(range),1));
        end
        
        function setAlign(self, aligns)
        % Set the alignment of each column
        %
        % Example:
        %   >>> T.setAlign('lcr')
            if length(aligns) ~= self.nColumns
                err('Error: aligns must have l,c, or r character for each column')
            end
            self.colAlign = aligns;
        end
        
        function setOmitRange(self, omitRange)
        % Set columns that can be omitted if there is not enough space
            self.omitRange = omitRange;
        end

        function addRow(self, row, loc, sep)
        % adds row to the table
        %
        % Convention: loc = 0 means add a row above the table and 
        %             loc = T.nrows adds a row below the table
        %
        % Example:
        %   >>> T = replab.str.Table(elements);
        %   >>> T.addRow({...}, 0)
            dim = size(self.elements);
            array = cell(dim(1) + 1, dim(2)); 
            array(loc + 1, :) = self.convertToChars(row);
            array(loc + 2:end, :) = self.elements(loc + 1:end, :);
            array(1:loc, :) = self.elements(1:loc, :);
            self.elements = array;
            if ~isempty(self.rowSep)
                self.rowSep = [self.rowSep(1:loc), {sep}, self.rowSep(loc+1:end)];
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
        % Example:
        %   >>> T = replab.str.Table(elements);
        %   >>> T.addColumn({...}, 0)
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
                                    columnwidth + length(sep), self.columnLengths(loc+1:end)];
        end

        function addColumnNames(self, colNames)
        % Adds row names to the current table
        % Convention: should not include an entry above column names
            if self.colName
                self.elements = self.elements(2:end,:);
                self.rowSep = self.rowSep(2:end);
            end
            if self.rowName
                self.addRow([{''}, colNames], 0);
            else
                self.addRow(colNames, 0);
            end
            self.colName = 1;
        end

        function addRowNames(self, rowNames)
        % Adds column names to the current table
        % Convention: should not include an entry before row names
            if self.rowName
                self.elements = self.elements(:, 2:end);
                self.colAlign = self.colAlign(2:end);
                self.colSep = self.colSep(2:end);
                self.columnLengths = self.columnLengths(2:end);
            end
            if self.colName
                self.addColumn([{''}, rowNames], 0, '  ', 'c');
            else
                self.addColumn(rowNames, 0, '  ', 'c');
            end
            self.rowName = 1;
        end


    end
    
    methods (Static)
        
        function chararray = convertToChars(elements)
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
            diff = omitRange(2:end) - omitRange(1:end-1);
            f = find(diff ~= 1);
            ellipsisLocs = omitRange([1, f + 1]);
            omitRange([1,f + 1]) = [];
            hiddenRange = omitRange;
        end
       
        
    end

end