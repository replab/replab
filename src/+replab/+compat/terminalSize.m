function [nRows nCols] = terminalSize
% Retrieves the size of the current terminal window
%
% Returns ``[]`` for both the number of rows and columns if that information is not available
%
% Returns
% -------
% nRows:
%   integer or ``[]``: Number of rows
% nCols:
%   integer or ``[]``: Number of columns
    nRows = [];
    nCols = [];
    if replab.compat.isOctave
        try
            pair = terminal_size;
            nRows = pair(1);
            nCols = pair(2);
        catch
        end
    else
        try
            pair = matlab.desktop.commandwindow.size;
            nRows = pair(2);
            nCols = pair(1);
        catch
        end
    end
    if isempty(nRows) || isempty(nCols)
        % to be sure
        nRows = [];
        nCols = [];
    end
end
