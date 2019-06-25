function [lines overLimit] = longStr(obj, maxRows, maxColumns)
% Returns a multiline description of the given object, that fits within the given width/height limit
% (this is the fallback implementation; user  'replab.longStr'
%
%        obj: Object to pretty print
%
%    maxRows: maximum number of rows
% maxColumns: maximum number of columns
%
% Returns a cell array of strings 's' and a Boolean 'overLimit' that states whether the output has been shortened
% because it ran over limit.
%   
% Thank you Matlab for being a regular language, we totally don't have to consider a thousand
% particular cases below.
    header = [];
    body = {};
    overLimit = false;
    if isscalar(obj) && (isobject(obj) || isstruct(obj))
        [names values] = replab.str.fieldsList(obj);
        n = length(names);
        header = [class(obj) ' instance with ' replab.str.pluralize(n, 'field') ':'];
        table = cell(n, 3);
        for i = 1:n
            name = names{i};
            value = values{i};
            table{i,1} = name;
            table{i,2} = ': ';
            table{i,3} = replab.shortStr(value, maxColumns - length(name) - 2);
        end
        body = replab.str.align(table, 'rcl');
    elseif isscalar(obj)
        header = replab.shortStr(obj, maxColumns);
        if iscell(obj)
            header = ['{' header '}'];
        end
        body = {};
    elseif isvector(obj)
        [header overShortLimit] = replab.shortStr(obj, maxColumns);
        if overShortLimit
            header = replab.str.tinyStr(obj);
            n = length(obj);
            if n >= maxRows
                n = maxRows - 1;
                overLimit = true;
            else
                overLimit = false;
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
            body{end+1, 1} = '...';
        end
    elseif ismatrix(obj)
        [header overShortLimit] = replab.shortStr(obj, maxColumns);
        if numel(obj) > 9 || overShortLimit
            header = replab.str.tinyStr(obj);
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
        warning('Pretty printing not implemented')
        header = replab.str.tinyStr(obj);
    end
    for i = 1:length(body)
        if length(body{i}) > maxColumns
            b = body{i};
            body{i} = [b(1:maxColumns-3) '...'];
            overLimit = true;
        end
    end
    if isequal(header, [])
        lines = body;
    else
        if length(header) > maxColumns
            overLimit = true;
            header = replab.str.tinyStr(obj);
        end
        lines = vertcat({header}, body);
    end
end
