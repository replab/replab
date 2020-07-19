function s = align(table, spec)
% Renders a table using the given LaTeX-inspired specification
%
% Args:
%   table (cell(\*,\*) of charstring): Cell contents as strings, cell content cannot be multiline
%   spec (char(1,\*)): Contains as many 'l', 'c', or 'r' to prescribe alignement as there are columns
%
% Returns:
%   (cell(\*,1) of charstring): Formatted text lines
    T = replab.str.Table(table, 'colSep', '');
    T.setAlign(1:T.nColumns, spec);
    s = T.format(1e6, 1e6);
    s = strsplit(s, '\n');
    s = s(:);
    s = cellfun(@deblank, s, 'uniform', 0);
end
