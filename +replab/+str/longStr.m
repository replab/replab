function lines = longStr(obj, maxRows, maxColumns)
% Returns a multiline description of the given object, that fits within the given width/height limit
% (this is the fallback implementation; call 'replab.longStr' in user code)
%
%        obj: Object to pretty print
%
%    maxRows: maximum number of rows
% maxColumns: maximum number of columns
%
% Returns a nRows x 1 cell array of strings 's'.
%   
% Thank you Matlab for being a regular language, we totally don't have to consider a thousand
% particular cases below.
    header = [];
    body = {};
    if isscalar(obj) && (isobject(obj) || isstruct(obj))
        [names values] = replab.str.fieldsList(obj);
        n = length(names);
        header = replab.headerStr(obj);
        table = cell(n, 3);
        for i = 1:n
            name = names{i};
            value = values{i};
            maxLength = maxColumns - length(name) - 2;
            str = replab.shortStr(value, maxLength);
            if length(str) > maxLength
                str = replab.headerStr(value);
            end
            table{i,1} = name;
            table{i,2} = ': ';
            table{i,3} = str;
        end
        body = replab.str.align(table, 'rcl');
    elseif isscalar(obj)
        header = replab.shortStr(obj, maxColumns);
        if iscell(obj)
            header = ['{' header '}'];
        end
        body = {};
    elseif isvector(obj)
        header = replab.shortStr(obj, maxColumns);
        if length(header) > maxColumns
            header = replab.str.headerStr(obj);
            n = length(obj);
            if n >= maxRows
                n = maxRows + 1; % do not print more, but keep extra stuff to signal overflow
            end
            table = cell(n, 3);
            for i = 1:n
                table{i,1} = ' ';
                table{i,2} = replab.str.cellStr(obj, maxColumns - 2, i);
                table{i,3} = ' ';
            end
            [lb rb] = replab.str.brackets(obj);
            table{1,1} = lb;
            table{n,3} = rb;
            body = replab.str.align(table, 'clc');
        end
    elseif ismatrix(obj)
        header = replab.shortStr(obj, maxColumns);
        if numel(obj) > 9 || length(header) > maxColumns
            header = replab.str.headerStr(obj);
            nR = size(obj, 1)
            nC = size(obj, 2)
            table = cell(nR, nC*2 + 1);
            for r = 1:nR
                table{r, 1} = ' ';
                for c = 1:nC
                    table{r, 2*c} = replab.str.cellStr(obj, maxColumns - 2, r, c);
                    table{r, 2*c+1} = ' ';
                end
            end
            [lb rb] = replab.str.brackets(obj);
            table{1, 1} = lb;
            table{nR, nC*2 + 1} = rb;
            body = replab.str.align(table, repmat('l', [1 nC*2+1]));
        end
    else
        header = replab.str.shortStr(obj);
    end
    if isequal(header, [])
        lines = body;
    else
        if length(header) > maxColumns
            header = replab.str.headerStr(obj);
        end
        lines = vertcat({header}, body);
    end
end
