function s = align(table, spec)
% Renders a table using the given LaTeX-inspired specification
%
% spec: contains as many 'l', 'c', or 'r' to prescribe alignement as there are columns
    nR = size(table, 1);
    nC = size(table, 2);
    s = cell(nR, 1);
    W = max(cellfun(@(x) length(x), table));
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
