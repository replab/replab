function [lines overLimit] = longFit(lines, overLimit, maxRows, maxColumns)
% Fits the string description given in 'lines' in the box of size maxRows x maxColumns
%
%      lines: column cell array where each cell element is a string representing a row
%  overLimit: if we had to trim due to breaking a limit
%             the return value is an 'or' of the previous 'overLimit' and whatever cut
%             is done here
%    maxRows: maximum number of rows
% maxColumns: maximum number of columns
    if length(lines) > maxRows
        lines = lines(1:maxRows-1);
        lines{maxRows} = '...';
        overLimit = true;
    end
    for i = 1:length(lines)
        if length(lines{i}) > maxColumns
            l = deblank(lines{i});
            if length(l) > maxColumns
                lines{i} = [l(1:maxColumns-3) '...'];
                overLimit = true;
            else
                lines{i} = [l repmat(' ', [1 maxColumns-length(l)])];
            end
        end
    end
end
