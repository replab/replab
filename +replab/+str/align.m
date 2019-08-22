function s = align(table, spec)
% Renders an aligned table using the given alignment of text in columns
%
% Args:
%   table (cell 2D array of strings): 2D table of data to pretty print
%   spec (char row vector): contains as many ``'l'``, ``'c'``, or ``'r'`` 
%                           to prescribe alignement as there are columns in the table
%
% Returns:
%   A column cell vector containing the lines of the table
    nR = size(table, 1);
    nC = size(table, 2);
    s = cell(nR, 1);
    Wtable = cellfun(@(x) length(x), table);
    W = max(Wtable, [], 1);
    for r = 1:nR
        row = char([]);
        for c = 1:nC
            el = table{r, c};
            pad = W(c) - length(el);
            switch spec(c)
              case 'l'
                lp = 0;
                rp = pad;
              case 'c'
                lp = ceil(pad/2);
                rp = pad - lp;
              case 'r'
                lp = pad;
                rp = 0;
            end
            row = [row repmat(' ', [1 lp]) el repmat(' ', [1 rp])];
        end
        s{r} = row;
    end
end
