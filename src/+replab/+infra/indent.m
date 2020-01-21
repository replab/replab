function indentedLines = indent(lines, nbSpaces, cellMode)
% indents lines
%
% Args:
%     lines (charstring or cell array of charstring): the text to be indented
%     nbSpaces (integer): the number of spaces to add at the beginning
%     cellMode (boolean, optional): whether the input should be expected to
%         be a cell array of charstring or a charstring. Default is true.
%
% Returns:
%     charstring or cell array of charstring: The indented text. The type
%         is a cell array iff parameter ``cellMode`` is true.
%
% Example:
%      >>> replab.infra.indent({'Some text'}, 2);
%        ans =
%        1x1 cell array
%        {'  Some txt'}

    if nargin < 3
        cellMode = true;
    end

    if iscell(lines) ~= cellMode
        if cellMode
            error('Lines should be provided in a cell array');
        else
            error('Lines should be provided as a charstring');
        end
    end

    if ~cellMode
        lines = {lines};
    end

    indentation = char(32*ones(1,nbSpaces));
    indentedLines = cell(size(lines));
    for i = 1:numel(lines)
        indentedLines{i} = [indentation, lines{i}];
        indentedLines{i} = strrep(indentedLines{i},char(10),[char(10), indentation]);
    end

    if ~cellMode
        indentedLines = indentedLines{1};
    end

end