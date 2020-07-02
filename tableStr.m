function tbstr = tableStr(Entries, RowNames, ColumnNames)
% Gives string that will display a pretty table
%
% Compatible with Octave
%
% Args:
%   Entries (cell array): body of the table
%   RowNames (optional cell array): top row of table
%   ColumnNames (optional cell array): left column of table
%
% Returns:
%   t (vector of strings): input to display to show table
    
    
    dim = size(Entries);
    char_arr = cell(dim);
    for i = 1:dim(1)
        for j = 1:dim(2)
            elmt = Entries{i, j};
            if ischar(elmt)
                charac = elmt;
            elseif isstring(elmt) % will only be true in MATLAB
                charac = convertStringsToChars(elmt);
            else
                charac = replab.shortStr(elmt);
            end
            char_arr{i, j} = charac;
        end
    end
            
%%% This worked for MATLAB
%     char_arr = cellfun(@change_to_char, Entries, 'UniformOutput',false);

    if nargin == 1
        rows = 0; cols = 0;
    elseif nargin == 2
        rows = 1; cols = 0;
    else
        rows = 1; cols = 1;
    end
    if rows
        char_arr(2:dim(1) + 1, :) = char_arr;
        for i = 1:dim(1)
            elmt = RowNames{i};
            if ischar(elmt)
                charac = elmt;
            elseif isstring(elmt) % will only be true in MATLAB
                charac = convertStringsToChars(elmt);
            else
                charac = replab.shortStr(elmt);
            end
            char_arr{1, i} = charac;
        end
        if cols
            char_arr(:, 2:dim(2) + 1) = char_arr;
            char_arr(1, 1) = {' '};
            for i = 1:dim(2)
                elmt = ColumnNames{i};
                if ischar(elmt)
                    charac = elmt;
                elseif isstring(elmt) % will only be true in MATLAB
                    charac = convertStringsToChars(elmt);
                else
                    charac = replab.shortStr(elmt);
                end
                char_arr{i + 1, 1} = charac;
            end
        end
    end

    max_len = max(max(cellfun(@length, char_arr))) + 2;
    dim = size(char_arr);
    max_term_len = floor(replab.settings.strMaxColumns / dim(2));
    max_len = min(max_len, max_term_len);
    
    for i = 1:dim(1)
        for j = 1:dim(2)
            elmt = char_arr{i, j};
            if length(elmt) > max_len - 2
                sized_elmt = [' ', elmt(1:max_len-3), '... '];
            else
                padding = ceil((max_len - 2 - length(elmt)) / 2);
                sized_elmt = [' ', repmat(' ', 1, max_len - 2 - length(elmt) - padding), elmt, ...
                              repmat(' ', 1, padding), ' '];
            end
            char_arr{i, j} = sized_elmt;
        end
    end
    
    char_arr(:, dim(2) + 1) = cell(1, dim(1));
    for i = 1:dim(1)
        char_arr{i, dim(2) + 1} = char(10);
    end
    tbstr = cell2mat(char_arr);

%%% Octave can't use nested functions
%     function char = change_to_char(elmt)
%         if ischar(elmt)
%             char = elmt;
%         elseif isstring(elmt) % will only be true in MATLAB
%             char = convertStringsToChars(elmt);
%         else
%             char = replab.shortStr(elmt);
%         end
%     end

end